function [hard_out,soft_out]=decoder_all_algorithm(in,position,num_iterate,algorithm,Lc)
%****************************************************************
% 内容概述：3种算法turbo解码器,in是RSC编码器输出
%          利用硬件化的方式实现TURBO码的译码
%          生成矩阵按照3GPP标准为[1 1 0 1;1 0 1 1]
%          输入为经过高斯信道的RSC软输入，而输出为软、硬输出
%****************************************************************
m=4;    %尾比特个数
g=[1 0 1 1;1 1 0 1];
L_seq=length(in);
in1=in(1:2,:);
in2=in(3:4,:);

e_p=zeros(1,L_seq);

for it=1:num_iterate
    a_p(position)=e_p(1:L_seq-m);  %解交织
    a_p(L_seq-m+1:L_seq)=0;  %尾比特部分不计算外部信息
    switch algorithm
        case 1
            [so,e_p] = constituent_decoder_logmap(in1,a_p,Lc); %log-map译码
        case 2
            [so,e_p] = constituent_decoder_max(in1,a_p,Lc); %max-log-map译码
        case 3
            [so,e_p] = constituent_decoder_map(in1,a_p,Lc); %map译码
        case 4
            [so,e_p] = constituent_decoder_sova(in1,g,a_p,1); 
    end
    
    a_p(1:L_seq-m)=e_p(position);  %交织
    a_p(L_seq-m+1:L_seq)=0;  %尾比特部分不计算外部信息
    switch algorithm
        case 1
            [so,e_p] = constituent_decoder_logmap(in2,a_p,Lc); %log-map译码
        case 2
            [so,e_p] = constituent_decoder_max(in2,a_p,Lc); %max-log-map译码
        case 3
            [so,e_p] = constituent_decoder_map(in1,a_p,Lc); %map译码
        case 4
            [so,e_p] = constituent_decoder_sova(in2,g,a_p,2); 
    end
end
% 解码结束，输出--------------------
soft_out(position)=so(1:L_seq-m); %输出解交织
hard_out=(sign(soft_out)+1)/2;  %取整为0、1形式
end