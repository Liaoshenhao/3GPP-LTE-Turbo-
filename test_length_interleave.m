%****************************************************************
% ���ݸ�����TURBO��������֮��֯����AWGN�ŵ����ԣ��˵�ʽ��
%          �ﵽ��֡�޼���ֹͣ��ǰSNR��Ĳ��ԣ���ʡ������
%          MAX-LOG-MAP�㷨��5����������20��֡��
%          ���鰴���ϲ������У��ڻ�����ƽ�����ߵ�ͬʱ��ʡʱ��
%          �˴�������ʱ���Ű��죨��������֡�ޣ�
%****************************************************************
clc;
clear all;

algorithm = input('�����㷨��1:LOG-MAP,2:MAX-LOG-MAP��Ĭ�ϣ�,3:MAP��ȱʡ��,4:SOVA��:');
if isempty(algorithm)
   algorithm =2;
end

length_interleave = input('��֯����=֡����β���س��ȡ�512 1024 1536��:');
if isempty(length_interleave)
   length_interleave = [512 1024 1536];
end

iteration = input('����������5��:');
if isempty(iteration)
   iteration = 5;
end

ferrlim = input('��֡��(�ﵽ���޼���ֹͣ��ǰSNR��Ĳ���)��20��:');
if isempty(ferrlim)
   ferrlim =20;
end

max_EbNo = input('�������Eb/No��dB����2.0��:');
if isempty(max_EbNo)
   max_EbNo =2.0;
end
 
step_EbNo = input('����Eb/No������dB����0.4��:');
if isempty(step_EbNo)
   step_EbNo=0.4;
end

save_mat = input('�Ƿ񱣴��������MAT�ļ� ��1-���棨ȱʡ��,0-�����桿:');
if isempty(save_mat)
   save_mat=1;
end

if save_mat==1
    matFileName = input('MAT�ļ��� ��''data2.mat''��:');
    if isempty(matFileName)
        matFileName='data2.mat';
    end
end


time_begin=datestr(now);
rate=1/3;           %����
m=4;                    %β������
fading_a=1;             %Fading amplitude
EbNo=0:step_EbNo:max_EbNo;                            %EbNo�Ĳ�����
EbNoLinear=10.^(EbNo.*0.1);     %transfer dB to linear

err_counter=zeros(length(length_interleave),length(EbNo));        %��ʼ��������ؼ�����
nferr= zeros(length(length_interleave),length(EbNo));             %��ʼ������֡������
ber=zeros(length(length_interleave),length(EbNo));                 %��ʼ�����������

for ii=1:length(length_interleave) %��֯����
    random_in=round(rand(1,length_interleave(ii)));  %�����
    [turbo_out,position]=encoder(random_in);      %����
    num_block_size=length_interleave(ii)+m;     %���ԵĿ�ߴ磬ָ����β���ص�������ϵͳϵ�г���

    for nEN=1:length(EbNo) %EbNo�Ĳ�����
        %L_c=4*fading_a*EbNoLinear(nEN)*rate; %�ŵ����Ŷ�
        L_c=1;
        sigma=1/sqrt(2*rate*EbNoLinear(nEN));
        nframe = 0;    % clear counter of transmitted frames
        
        
        if nEN==1 || ber(ii,nEN-1)>9.0e-6
            while nferr(ii,nEN)<ferrlim        %nferr:��ǰ����������EbNo��Ĵ���֡��
                nframe = nframe + 1; 
                noise=randn(4,num_block_size);    %����
                soft_in=L_c*(turbo_out+sigma*noise);     %��Ϣ��������
                %soft_in=awgn(turbo_out,15);     %awgn noise
                [hard_out,soft_out]=decoder_all_algorithm(soft_in,position,iteration,algorithm,L_c); %����
                errs=length(find(hard_out(1:length_interleave(ii))~=random_in));%��ǰ�����bit��
                
                if errs>0 
                    err_counter(ii,nEN)=err_counter(ii,nEN)+errs;
                    nferr(ii,nEN)=nferr(ii,nEN)+1;
                end
                fprintf('��ǰEbNo�㣺%1.2fdB���Ѽ��㣺%2.0f֡�����У�%2.0f��֡\n',...
                    EbNo(nEN),nframe,nferr(ii,nEN));
            end
            ber(ii,nEN) = err_counter(ii,nEN)/nframe/(length_interleave(ii));%�������
        else
            ber(ii,nEN)=NaN;
        end
        fprintf('��֯���ȣ�%1.0f��EbNo��%1.2fdB�������ʣ�%8.4e��\n',...
            length_interleave(ii),EbNo(nEN),ber(ii,nEN));
        if save_mat==1
            save (matFileName,'EbNo','ber');
        end
    end
end

semilogy(EbNo,ber(1,:),EbNo,ber(2,:),EbNo,ber(3,:));
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
title('3GPP��׼ Max-Log-MAP�����㷨,5��������,1/3����');
legend('512bit','1024bit','1536bit');

time_end=datestr(now);

fprintf('------------------��ϲ�㣡������ɣ�--------------------\n'); 
disp([' ������ʼʱ��:',time_begin,'=>',time_end])
fprintf(' ��������=%2dbit����֯����=%4d %4d %4d\n',iteration,length_interleave(1),length_interleave(2),length_interleave(3));

fprintf(' �������Eb/No = %2.1fdB������Eb/No���� = %2.1fdB\n',max_EbNo,step_EbNo);
switch algorithm
    case 1
        fprintf(' �����㷨��LOG-MAP\n');
    case 2
        fprintf(' �����㷨��MAX-LOG-MAP\n');
    case 3
        fprintf(' �����㷨������MAX-LOG-MAP\n');
    case 4
        fprintf(' �����㷨��SOVA\n');
end
if save_mat==1
    fprintf(' ����������� = %4s\n',matFileName);
end    
fprintf('-------------------------------------------------------\n'); 