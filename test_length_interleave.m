%****************************************************************
% 内容概述：TURBO译码性能之交织长度AWGN信道测试（菜单式）
%          达到误帧限即可停止当前SNR点的测试，节省计算量
%          MAX-LOG-MAP算法，5迭代次数，20误帧限
%          建议按以上参数运行，在画出较平滑曲线的同时节省时间
%          此代码运行时间大概半天（可增大误帧限）
%****************************************************************
clc;
clear all;

algorithm = input('译码算法【1:LOG-MAP,2:MAX-LOG-MAP（默认）,3:MAP（缺省）,4:SOVA】:');
if isempty(algorithm)
   algorithm =2;
end

length_interleave = input('交织长度=帧长－尾比特长度【512 1024 1536】:');
if isempty(length_interleave)
   length_interleave = [512 1024 1536];
end

iteration = input('迭代次数【5】:');
if isempty(iteration)
   iteration = 5;
end

ferrlim = input('误帧限(达到此限即可停止当前SNR点的测试)【20】:');
if isempty(ferrlim)
   ferrlim =20;
end

max_EbNo = input('测试最大Eb/No（dB）【2.0】:');
if isempty(max_EbNo)
   max_EbNo =2.0;
end
 
step_EbNo = input('测试Eb/No步长（dB）【0.4】:');
if isempty(step_EbNo)
   step_EbNo=0.4;
end

save_mat = input('是否保存仿真结果到MAT文件 【1-保存（缺省）,0-不保存】:');
if isempty(save_mat)
   save_mat=1;
end

if save_mat==1
    matFileName = input('MAT文件名 【''data2.mat''】:');
    if isempty(matFileName)
        matFileName='data2.mat';
    end
end


time_begin=datestr(now);
rate=1/3;           %码率
m=4;                    %尾比特数
fading_a=1;             %Fading amplitude
EbNo=0:step_EbNo:max_EbNo;                            %EbNo的采样点
EbNoLinear=10.^(EbNo.*0.1);     %transfer dB to linear

err_counter=zeros(length(length_interleave),length(EbNo));        %初始化错误比特计数器
nferr= zeros(length(length_interleave),length(EbNo));             %初始化错误帧计数器
ber=zeros(length(length_interleave),length(EbNo));                 %初始化错误比特率

for ii=1:length(length_interleave) %交织长度
    random_in=round(rand(1,length_interleave(ii)));  %随机数
    [turbo_out,position]=encoder(random_in);      %编码
    num_block_size=length_interleave(ii)+m;     %测试的块尺寸，指包含尾比特的软输入系统系列长度

    for nEN=1:length(EbNo) %EbNo的采样点
        %L_c=4*fading_a*EbNoLinear(nEN)*rate; %信道置信度
        L_c=1;
        sigma=1/sqrt(2*rate*EbNoLinear(nEN));
        nframe = 0;    % clear counter of transmitted frames
        
        
        if nEN==1 || ber(ii,nEN-1)>9.0e-6
            while nferr(ii,nEN)<ferrlim        %nferr:当前迭代次数、EbNo点的错误帧数
                nframe = nframe + 1; 
                noise=randn(4,num_block_size);    %噪声
                soft_in=L_c*(turbo_out+sigma*noise);     %信息噪声叠加
                %soft_in=awgn(turbo_out,15);     %awgn noise
                [hard_out,soft_out]=decoder_all_algorithm(soft_in,position,iteration,algorithm,L_c); %译码
                errs=length(find(hard_out(1:length_interleave(ii))~=random_in));%当前点错误bit数
                
                if errs>0 
                    err_counter(ii,nEN)=err_counter(ii,nEN)+errs;
                    nferr(ii,nEN)=nferr(ii,nEN)+1;
                end
                fprintf('当前EbNo点：%1.2fdB；已计算：%2.0f帧；其中：%2.0f误帧\n',...
                    EbNo(nEN),nframe,nferr(ii,nEN));
            end
            ber(ii,nEN) = err_counter(ii,nEN)/nframe/(length_interleave(ii));%误比特率
        else
            ber(ii,nEN)=NaN;
        end
        fprintf('交织长度：%1.0f；EbNo：%1.2fdB；误码率：%8.4e；\n',...
            length_interleave(ii),EbNo(nEN),ber(ii,nEN));
        if save_mat==1
            save (matFileName,'EbNo','ber');
        end
    end
end

semilogy(EbNo,ber(1,:),EbNo,ber(2,:),EbNo,ber(3,:));
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
title('3GPP标准 Max-Log-MAP译码算法,5迭代次数,1/3码率');
legend('512bit','1024bit','1536bit');

time_end=datestr(now);

fprintf('------------------恭喜你！测试完成！--------------------\n'); 
disp([' 仿真起始时间:',time_begin,'=>',time_end])
fprintf(' 迭代次数=%2dbit；交织长度=%4d %4d %4d\n',iteration,length_interleave(1),length_interleave(2),length_interleave(3));

fprintf(' 测试最大Eb/No = %2.1fdB；测试Eb/No步长 = %2.1fdB\n',max_EbNo,step_EbNo);
switch algorithm
    case 1
        fprintf(' 译码算法：LOG-MAP\n');
    case 2
        fprintf(' 译码算法：MAX-LOG-MAP\n');
    case 3
        fprintf(' 译码算法：门限MAX-LOG-MAP\n');
    case 4
        fprintf(' 译码算法：SOVA\n');
end
if save_mat==1
    fprintf(' 保存仿真结果到 = %4s\n',matFileName);
end    
fprintf('-------------------------------------------------------\n'); 