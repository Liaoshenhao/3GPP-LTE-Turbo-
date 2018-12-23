function [turbo_out,position]=encoder(in)
%****************************************************************
% 内容概述：Turbo码编码器。
%          生成矩阵按照3GPP标准为[1 1 0 1;1 0 1 1]
%          输入为随机生成的系统信息位，而输出为BPSK调制信号和交织位置
%****************************************************************

k=length(in);
Ck1=in;
x1 = zeros(1,k+4);x2 = zeros(1,k+4);
z1 = zeros(1,k+4);z2 = zeros(1,k+4);

%first RSC encoder
reg1=zeros(1,3); %初始化移位寄存器
for i=1:k
    a = mod((Ck1(1,i)+reg1(1,2)+reg1(1,3)),2); %前K bit移位寄存器的更新值
    x1(1,i) = Ck1(1,i);
    z1(1,i) = mod((a+reg1(1,1)+reg1(1,3)),2);
    reg1 = [a,reg1(1,[1:2])];
end
%trellis termination for RSC1
for m=1:3
    x1(1,k+m) = mod((reg1(1,2)+reg1(1,3)),2);
    z1(1,k+m) = mod((reg1(1,1)+reg1(1,3)),2);
    reg1 = [0,reg1(1,[1:2])];      %尾比特用0更新寄存器
end

%交织
position=interleaver_3GPP(k);
Ck2=Ck1(position);
x2(1,1:k)=Ck2;

%second RSC encoder
reg2=zeros(1,3); %initial register
for j=1:k
    d = mod((Ck2(1,j)+reg2(1,2)+reg2(1,3)),2);
    z2(1,j) = mod((d+reg2(1,1)+reg2(1,3)),2);
    reg2 = [d,reg2(1,[1:2])];
end
%trellis termination for RSC2
for m=1:3
    x2(1,k+m) = mod((reg2(1,2)+reg2(1,3)),2);
    z2(1,k+m) = mod((reg2(1,1)+reg2(1,3)),2);
    reg2 = [0,reg2(1,[1:2])];
end

%trellis termination for turbor encoder
%网格终止，参照标准
d0(1,1:k) = x1(1,1:k);
d1(1,1:k) = z1(1,1:k);
d2(1,1:k) = z2(1,1:k);
d3(1,1:k) = Ck2; 
d0(1,k+1)=x1(1,k+1);d0(1,k+2)=z1(1,k+2);d0(1,k+3)=x2(1,k+1);d0(1,k+4)=z2(1,k+2);
d1(1,k+1)=z1(1,k+1);d1(1,k+2)=x1(1,k+3);d1(1,k+3)=z2(1,k+1);d1(1,k+4)=x2(1,k+3);
d2(1,k+1)=x1(1,k+2);d2(1,k+2)=z1(1,k+3);d2(1,k+3)=x2(1,k+2);d2(1,k+4)=z2(1,k+3);
d3(1,k+1:k+4)=d0(1,k+1:k+4);%尾比特不交织
turbo_out(1,:)=d0;% 信息位
turbo_out(2,:)=d1;% 校验位
turbo_out(3,:)=d3;% 交织后信息位 
turbo_out(4,:)=d2;% 交织后校验位
turbo_out=2*turbo_out-ones(size(turbo_out));
% 调制 将 1 调制成  +1
%         0         -1



