% ����


%%                             
%% ����PRN ��ȡ��Ƶ����
clear
[codeai, codeaq, codebi, codebq] = genPRN(7);
%
clear cor_sum_ai cor_sum_aq cor_sum_bi cor_sum_bq
close all
ifdata = ReadIF('./20210428_A40ms.txt');
ifdata = ifdata(22581:end);
%% fft����
clear cor_sum_ai cor_sum_aq cor_sum_bi cor_sum_bq
% [cor_sum_ai, doppler_ai, codep_ai] = acquire_fft(ifdata, codeai, 1);                       
% [cor_sum_aq, doppler_aq, codep_aq] = acquire_fft(ifdata, codeaq, 1);
% [cor_sum_bi, doppler_bi, codep_bi] = acquire_fft(ifdata, codebi, 0);
[cor_sum_bq, doppler_bq, codep_bq] = acquire_fft(ifdata, codebq, 0);
%% ��ȷ����
clear cor_sum_ai cor_sum_aq cor_sum_bi cor_sum_bq
freq_acq = 920;
code_acq = 1;
% [cor_sum_ai] = acquire_acc(ifdata, codeai, 1, freq_acq, code_acq);
% [cor_sum_aq] = acquire_acc(ifdata, codeaq, 1, freq_acq, code_acq);
% [cor_sum_bi] = acquire_acc(ifdata, codebi, 0, freq_acq, code_acq);
[cor_sum_bq] = acquire_acc(ifdata, codebq, 0, freq_acq, code_acq);
%%
fun1 = @(x) surf(((1:25)*400-5000), ((1:58000)/58000*10230), (x(1:20:end,1:2:end).*conj(x(1:20:end,1:2:end))).','LineStyle', 'none');
%%
figure(1)
fun1(cor_sum_ai)
xlabel('freq (Hz)')
ylabel('code')
%%
figure(2)
fun1(cor_sum_aq)
xlabel('freq (Hz)')
ylabel('code')
%%
figure(3)
fun1(cor_sum_bi)
xlabel('freq (Hz)')
ylabel('code')
%%
fun1 = @(x) surf(((1:25)*400-5000), ((1:58000)/58000*10230), (x(1:20:end,1:2:end).*conj(x(1:20:end,1:2:end))).','LineStyle', 'none');
figure(4)
fun1(cor_sum_bq)
xlabel('freq (Hz)')
ylabel('code')
%%
acqui_cp = 93421;
fun1 = @(x) surf( (x(1:1:end,acqui_cp-20:acqui_cp+20).*conj(x(1:1:end,acqui_cp-20:acqui_cp+20))).','LineStyle', 'none');
figure(1)
fun1(cor_sum_ai)
xlabel('freq (Hz)')
ylabel('code')
figure(2)
fun1(cor_sum_aq)
xlabel('freq (Hz)')
ylabel('code')
figure(3)
fun1(cor_sum_bi)
xlabel('freq (Hz)')
ylabel('code')
figure(4)
fun1(cor_sum_bq)
xlabel('freq (Hz)')
ylabel('code')
%%
fun1 = @(x) surfc( (x.*conj(x)).','LineStyle', 'none');
figure(1)
fun1(cor_sum_ai)
xlabel('freq (Hz)')
ylabel('code')
figure(2)
fun1(cor_sum_aq)
xlabel('freq (Hz)')
ylabel('code')
figure(3)
fun1(cor_sum_bi)
xlabel('freq (Hz)')
ylabel('code')
figure(4)
fun1(cor_sum_bq)
xlabel('freq (Hz)')
ylabel('code')
%%
surf(cor_sum_aq(1:20:end,1:2:end).*conj(cor_sum_aq(1:20:end,1:2:end)),'LineStyle', 'none')
%% 
cor_pow = cor_sum.*conj(cor_sum);
image(cor_pow.', 'CDataMapping','scaled');
colorbar
[r, w] = find(cor_pow == max(cor_pow, [], 'all'));
%% 
figure(2)
cor_pow = cor_sum_bq.*conj(cor_sum_bq);
image(cor_pow.', 'CDataMapping','scaled');
colorbar
[r, w] = find(cor_pow == max(cor_pow, [], 'all'));
%%
surf(cor_sum_bq(1:20:end,1:2:end).*conj(cor_sum_bq(1:20:end,1:2:end)),'LineStyle', 'none')
%%
save('./data/cor_sum.mat', 'cor_sum');
%%
function [cor_sum, freq_acqui, code_phase] = acquire_fft(ifdata, code, lowband)
% @brief: ����������ؽ��
% @param:
% ifdata - ��Ƶ���� nx1 complex double
% code - Q·��Ƶ������
% lowband - lowband = 1 E5a�źţ�lowband = 0 E5b�ź�
% @return:
% cor_sum - �����ͬƵ���µ���ؽ��
% @warning

SystemClockFre = 5.8e7; %��¼�ط��ǲ�����58MHz
GALILEO_FREQ_COMPENSATE = -0.56439e6; %٤���Ե���ƵƵ��
GALILEO_CACODE_FREQ = 10.23e6; %٤���Ե�CA������10.23MHz
GALILEO_CACODE_NUM = 10230; %٤���Ե�CA�볤��

Freq_Range = 10e3;% �����Ƶ�ʷ�Χ20kHz -10kHz~10kHz
f0 = GALILEO_FREQ_COMPENSATE;% ���Ƶ���ƵƵ��
f_delta = 20;% ����Ƶ�ʲ����� 10Hz
fs = 1.5*GALILEO_CACODE_FREQ;% ���Ƶĸ��ز�Ƶ
cor_1ms = SystemClockFre/1000;% 1ms�ĵ��� 58000

cor_sum = zeros(Freq_Range/f_delta, cor_1ms*2);% ��ؽ�������С
index = ceil(((1:cor_1ms)/cor_1ms*GALILEO_CACODE_NUM));% ��prn���10230��ת����58000����ת������
code_1ms = code(index);
code_1ms = [code_1ms.' (zeros(cor_1ms,1)).'].';% ����1ms
code_1ms_F = fft(code_1ms);% fft
% �ж��Ƿ����±ߴ�
if lowband == 1
    f_IF = f0 - fs;
else
    f_IF = f0 + fs;
end
ifdata_2ms = ifdata(1:cor_1ms*2);% �ضϵ�2ms
% fft��ѭ�����
parfor i = 1:Freq_Range/f_delta
    f_temp = f_IF - Freq_Range/2 + i*f_delta;% e5aƵ�ʹ���ֵ
    DDCData = ifdata_2ms.*exp(-1i*2*pi*f_temp/SystemClockFre*(1:cor_1ms*2)).';%  �±�Ƶ����Ƶ
    DDCData_F = fft(DDCData);
    cor_sum(i,:) = (ifft(DDCData_F.*conj(code_1ms_F)));
end
% �󲶻��Ƶ�ʺ�����λ
cor_pow = cor_sum.*conj(cor_sum);
[r, w] = find(cor_pow == max(cor_pow, [], 'all'));
freq_acqui = - Freq_Range/2 + r*f_delta;% ����Ķ�����Ƶ��
code_phase = w/cor_1ms*GALILEO_CACODE_NUM;% ����λ

end


function [cor_sum] = acquire_acc(ifdata, code, lowband, freq_acq, code_acq)
% @brief: ��ȷ����
% @param:
% ifdata - ��Ƶ���� nx1 complex double
% code - Q·��Ƶ������
% lowband - lowband = 1 E5a�źţ�lowband = 0 E5b�ź�
% @return:
% cor_sum - �����ͬƵ���µ���ؽ��
% @warning
SystemClockFre = 5.8e7; %��¼�ط��ǲ�����58MHz
GALILEO_FREQ_COMPENSATE = -0.56439e6; %٤���Ե���ƵƵ��
GALILEO_CACODE_FREQ = 10.23e6; %٤���Ե�CA������10.23MHz
GALILEO_CACODE_NUM = 10230; %٤���Ե�CA�볤��

Freq_Range = 100;% �����Ƶ�ʷ�Χ50Hz -25Hz~25Hz
Code_Range = 4;% �������Ƭ��Χ -2~+2
f0 = GALILEO_FREQ_COMPENSATE;% ���Ƶ���ƵƵ��
f_delta = 1;% ����Ƶ�ʲ����� 1Hz
code_delta = 0.1;
fs = 1.5*GALILEO_CACODE_FREQ;% ���Ƶĸ��ز�Ƶ��
cor_1ms = SystemClockFre/1000;% 1ms�ĵ��� 58000

cor_sum = zeros(Freq_Range/f_delta, Code_Range/code_delta);% ��ؽ�������С
cor_sum1 = zeros(Freq_Range/f_delta, 1);

% index = ceil(((1:cor_1ms)/cor_1ms*GALILEO_CACODE_NUM));% ��prn���10230��ת����58000����ת������
% code_1ms = code(index);
% code_1ms = [code_1ms.' (zeros(cor_1ms,1)-1).'].';% ����1ms
% code_1ms_F = fft(code_1ms);% fft
code_2ms = [code.' (zeros(GALILEO_CACODE_NUM, 1)).'].';% ����� 2ms��α�����

% �ж��Ƿ����±ߴ�
if lowband == 1
    f_IF = f0 - fs;
else
    f_IF = f0 + fs;
end
ifdata_2ms = ifdata(1:cor_1ms*2);% �ضϵ�2ms
% ��ѭ�����
for i = 1:Freq_Range/f_delta
    f_temp = f_IF - Freq_Range/2 + i*f_delta + freq_acq;% e5aƵ�ʹ���ֵ
    DDCData = ifdata_2ms.*exp(1i*2*pi*f_temp/SystemClockFre*(1:cor_1ms*2)).';%  �±�Ƶ����Ƶ
    cor_sum2 = zeros(Code_Range/code_delta, 1);
    for j = 1:Code_Range/code_delta
        index = ceil((((1:cor_1ms*2)+j*code_delta-Code_Range/2 - code_acq + 1)/cor_1ms*GALILEO_CACODE_NUM));% ��prn���20460��ת����58000����ת������
        index = mod(index -1, GALILEO_CACODE_NUM*2)+1;
        cor_sum2(j) = sum(DDCData.*code_2ms(index));
    end
    cor_sum(i, :) = cor_sum2.';
end
% �󲶻��Ƶ�ʺ�����λ
cor_pow = cor_sum.*conj(cor_sum);
[r, w] = find(cor_pow == max(cor_pow, [], 'all'));
freq_acqui = - Freq_Range/2 + r*f_delta;% ����Ķ�����Ƶ��
code_phase = w/cor_1ms*GALILEO_CACODE_NUM;% ����λ
end