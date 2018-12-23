function [L_all,L_e] = constituent_decoder_sova(in, g, L_a, ind_dec) 
%****************************************************************
% 内容概述：子译码器。
%          利用硬件化的方式实现TURBO码的SOVA译码
%          生成矩阵按照3GPP标准为[1 1 0 1;1 0 1 1]
%          输入为经过高斯信道的RSC软输入，而输出为软输出和外部信息
% Output:
%       L_all: log ( P(x=1|y) ) / ( P(x=-1|y) )
%       注：根据网上文档整改成一个脚本，此算法运行慢
%****************************************************************
L_total = length(L_a);
[n,K] = size(g); 
m = K - 1;
nstates = 2^m;
Infty = 1e10;
for i=1:L_total
    rec(1,2*(i-1)+1)=in(1,i);
    rec(2*i)=in(2,i);
end
rec_s=0.5*rec;

% SOVA window size. Make decision after 'delta' delay. Decide bit k when received bits
% for bit (k+delta) are processed. Trace back from (k+delta) to k. 
delta = 30;    

% Set up the trellis defined by g.

for state=1:nstates
   state1=state-1;
   for j = 1:length( state1 )
        for i = m:-1:1
            state_o(j,m-i+1) = fix( state1(j)/ (2^(i-1)) );
            state1(j) = state1(j) - state_o(j,m-i+1)*2^(i-1);
        end
   end
    state_vector = state_o;
   
   % when receive a 0
   d_k = 0;
   a_k = rem( g(1,:)*[0 state_vector]', 2 ); %a_k为1 X 4矩阵
   
   for i=1:n
        out_0(i) = g(i,1)*a_k;
        for j = 2:K
             out_0(i) = xor(out_0(i),g(i,j)*state_vector(j-1));
        end;
   end

   state_0 = [a_k, state_vector(1:m-1)];
   out_0(1) = 0;
  
   % when receive a 1
   d_k = 1;
   a_k = rem( g(1,:)*[1 state_vector]', 2 );
   
   for i=1:n
        out_1(i) = g(i,1)*a_k;
        for j = 2:K
             out_1(i) = xor(out_1(i),g(i,j)*state_vector(j-1));
        end;
   end

   state_1 = [a_k, state_vector(1:m-1)];
   out_1(1) = 1;
   next_out(state,:) = 2*[out_0 out_1]-1;
   
   [dummy1, m1] = size( state_0 );

    for i = 1:m1
         vect_0(i) = 2^(m1-i);
    end
    int_state_0 = state_0*vect_0';
    
    [dummy2, m2] = size( state_1 );

    for i = 1:m2
         vect_1(i) = 2^(m2-i);
    end
    int_state_1 = state_1*vect_1';
   next_state(state,:) = [(int_state_0+1) (int_state_1+1)];
end

% find out which two previous states can come to present state
last_state = zeros(nstates,2);
for bit=0:1
   for state=1:nstates
      last_state(next_state(state,bit+1), bit+1)=state;
      last_out(next_state(state, bit+1), bit*2+1:bit*2+2) ...
         = next_out(state, bit*2+1:bit*2+2);
   end 
end



% Initialize path metrics to -Infty
for t=1:L_total+1
   for state=1:nstates
      path_metric(state,t) = -Infty;
   end
end

% Trace forward to compute all the path metrics
path_metric(1,1) = 0;
for t=1:L_total
   y = rec_s(2*t-1:2*t);
   for state=1:nstates
      sym0 = last_out(state,1:2);
      sym1 = last_out(state,3:4);
      state0 = last_state(state,1);
      state1 = last_state(state,2);
      Mk0 = y*sym0' - L_a(t)/2 + path_metric(state0,t);
      Mk1 = y*sym1' + L_a(t)/2 + path_metric(state1,t);
      
      if Mk0>Mk1
         path_metric(state,t+1)=Mk0;
         Mdiff(state,t+1) = Mk0 - Mk1;
         prev_bit(state, t+1) = 0;
      else
         path_metric(state,t+1)=Mk1;
         Mdiff(state,t+1) = Mk1 - Mk0;
         prev_bit(state,t+1) = 1;
      end

   end
end
      
% For decoder 1, trace back from all zero state, 
% for decoder two, trace back from the most likely state
if ind_dec == 1
   mlstate(L_total+1) = 1;
else
   mlstate(L_total+1) = find( path_metric(:,L_total+1)==max(path_metric(:,L_total+1)) );
end

% Trace back to get the estimated bits, and the most likely path
for t=L_total:-1:1
   est(t) = prev_bit(mlstate(t+1),t+1);
   mlstate(t) = last_state(mlstate(t+1), est(t)+1);
end

% Find the minimum delta that corresponds to a compitition path with different info. bit estimation.       
% Give the soft output
for t=1:L_total
   llr = Infty;
   for i=0:delta
      if t+i<L_total+1
         bit = 1-est(t+i);
         temp_state = last_state(mlstate(t+i+1), bit+1);
         for j=i-1:-1:0
            bit = prev_bit(temp_state,t+j+1);
            temp_state = last_state(temp_state, bit+1);
         end
         if bit~=est(t) 
            llr = min( llr,Mdiff(mlstate(t+i+1), t+i+1) );
         end
      end
   end
   L_all(t) = (2*est(t) - 1) * llr;
end

L_e = L_all - 2*rec_s(1,1:2:2*L_total) - L_a;  % extrinsic info.