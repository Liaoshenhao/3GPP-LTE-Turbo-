function [hard_out,soft_out]=decoder_all_algorithm(in,position,num_iterate,algorithm,Lc)
%****************************************************************
% ���ݸ�����3���㷨turbo������,in��RSC���������
%          ����Ӳ�����ķ�ʽʵ��TURBO�������
%          ���ɾ�����3GPP��׼Ϊ[1 1 0 1;1 0 1 1]
%          ����Ϊ������˹�ŵ���RSC�����룬�����Ϊ��Ӳ���
%****************************************************************
m=4;    %β���ظ���
g=[1 0 1 1;1 1 0 1];
L_seq=length(in);
in1=in(1:2,:);
in2=in(3:4,:);

e_p=zeros(1,L_seq);

for it=1:num_iterate
    a_p(position)=e_p(1:L_seq-m);  %�⽻֯
    a_p(L_seq-m+1:L_seq)=0;  %β���ز��ֲ������ⲿ��Ϣ
    switch algorithm
        case 1
            [so,e_p] = constituent_decoder_logmap(in1,a_p,Lc); %log-map����
        case 2
            [so,e_p] = constituent_decoder_max(in1,a_p,Lc); %max-log-map����
        case 3
            [so,e_p] = constituent_decoder_map(in1,a_p,Lc); %map����
        case 4
            [so,e_p] = constituent_decoder_sova(in1,g,a_p,1); 
    end
    
    a_p(1:L_seq-m)=e_p(position);  %��֯
    a_p(L_seq-m+1:L_seq)=0;  %β���ز��ֲ������ⲿ��Ϣ
    switch algorithm
        case 1
            [so,e_p] = constituent_decoder_logmap(in2,a_p,Lc); %log-map����
        case 2
            [so,e_p] = constituent_decoder_max(in2,a_p,Lc); %max-log-map����
        case 3
            [so,e_p] = constituent_decoder_map(in1,a_p,Lc); %map����
        case 4
            [so,e_p] = constituent_decoder_sova(in2,g,a_p,2); 
    end
end
% ������������--------------------
soft_out(position)=so(1:L_seq-m); %����⽻֯
hard_out=(sign(soft_out)+1)/2;  %ȡ��Ϊ0��1��ʽ
end