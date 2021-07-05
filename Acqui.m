% 捕获


%%                             
%% 生成PRN 读取中频数据
clear
[codeai, codeaq, codebi, codebq] = genPRN(7);
%
clear cor_sum_ai cor_sum_aq cor_sum_bi cor_sum_bq
close all
ifdata = ReadIF('./20210428_A40ms.txt');
ifdata = ifdata(22581:end);
%% fft捕获
clear cor_sum_ai cor_sum_aq cor_sum_bi cor_sum_bq
% [cor_sum_ai, doppler_ai, codep_ai] = acquire_fft(ifdata, codeai, 1);                       
% [cor_sum_aq, doppler_aq, codep_aq] = acquire_fft(ifdata, codeaq, 1);
% [cor_sum_bi, doppler_bi, codep_bi] = acquire_fft(ifdata, codebi, 0);
[cor_sum_bq, doppler_bq, codep_bq] = acquire_fft(ifdata, codebq, 0);
%% 精确捕获
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
% @brief: 输出捕获的相关结果
% @param:
% ifdata - 中频数据 nx1 complex double
% code - Q路导频数据码
% lowband - lowband = 1 E5a信号；lowband = 0 E5b信号
% @return:
% cor_sum - 输出不同频率下的相关结果
% @warning

SystemClockFre = 5.8e7; %记录回放仪采样率58MHz
GALILEO_FREQ_COMPENSATE = -0.56439e6; %伽利略的中频频率
GALILEO_CACODE_FREQ = 10.23e6; %伽利略的CA码速率10.23MHz
GALILEO_CACODE_NUM = 10230; %伽利略的CA码长度

Freq_Range = 10e3;% 捕获的频率范围20kHz -10kHz~10kHz
f0 = GALILEO_FREQ_COMPENSATE;% 估计的中频频率
f_delta = 20;% 捕获频率步进量 10Hz
fs = 1.5*GALILEO_CACODE_FREQ;% 估计的副载波频
cor_1ms = SystemClockFre/1000;% 1ms的点数 58000

cor_sum = zeros(Freq_Range/f_delta, cor_1ms*2);% 相关结果分配大小
index = ceil(((1:cor_1ms)/cor_1ms*GALILEO_CACODE_NUM));% 将prn码从10230个转换到58000个的转换序列
code_1ms = code(index);
code_1ms = [code_1ms.' (zeros(cor_1ms,1)).'].';% 补零1ms
code_1ms_F = fft(code_1ms);% fft
% 判断是否是下边带
if lowband == 1
    f_IF = f0 - fs;
else
    f_IF = f0 + fs;
end
ifdata_2ms = ifdata(1:cor_1ms*2);% 截断到2ms
% fft求循环相关
parfor i = 1:Freq_Range/f_delta
    f_temp = f_IF - Freq_Range/2 + i*f_delta;% e5a频率估计值
    DDCData = ifdata_2ms.*exp(-1i*2*pi*f_temp/SystemClockFre*(1:cor_1ms*2)).';%  下变频到中频
    DDCData_F = fft(DDCData);
    cor_sum(i,:) = (ifft(DDCData_F.*conj(code_1ms_F)));
end
% 求捕获的频率和码相位
cor_pow = cor_sum.*conj(cor_sum);
[r, w] = find(cor_pow == max(cor_pow, [], 'all'));
freq_acqui = - Freq_Range/2 + r*f_delta;% 捕获的多普勒频率
code_phase = w/cor_1ms*GALILEO_CACODE_NUM;% 码相位

end


function [cor_sum] = acquire_acc(ifdata, code, lowband, freq_acq, code_acq)
% @brief: 精确捕获
% @param:
% ifdata - 中频数据 nx1 complex double
% code - Q路导频数据码
% lowband - lowband = 1 E5a信号；lowband = 0 E5b信号
% @return:
% cor_sum - 输出不同频率下的相关结果
% @warning
SystemClockFre = 5.8e7; %记录回放仪采样率58MHz
GALILEO_FREQ_COMPENSATE = -0.56439e6; %伽利略的中频频率
GALILEO_CACODE_FREQ = 10.23e6; %伽利略的CA码速率10.23MHz
GALILEO_CACODE_NUM = 10230; %伽利略的CA码长度

Freq_Range = 100;% 捕获的频率范围50Hz -25Hz~25Hz
Code_Range = 4;% 捕获的码片范围 -2~+2
f0 = GALILEO_FREQ_COMPENSATE;% 估计的中频频率
f_delta = 1;% 捕获频率步进量 1Hz
code_delta = 0.1;
fs = 1.5*GALILEO_CACODE_FREQ;% 估计的副载波频率
cor_1ms = SystemClockFre/1000;% 1ms的点数 58000

cor_sum = zeros(Freq_Range/f_delta, Code_Range/code_delta);% 相关结果分配大小
cor_sum1 = zeros(Freq_Range/f_delta, 1);

% index = ceil(((1:cor_1ms)/cor_1ms*GALILEO_CACODE_NUM));% 将prn码从10230个转换到58000个的转换序列
% code_1ms = code(index);
% code_1ms = [code_1ms.' (zeros(cor_1ms,1)-1).'].';% 补零1ms
% code_1ms_F = fft(code_1ms);% fft
code_2ms = [code.' (zeros(GALILEO_CACODE_NUM, 1)).'].';% 补零后 2ms的伪随机码

% 判断是否是下边带
if lowband == 1
    f_IF = f0 - fs;
else
    f_IF = f0 + fs;
end
ifdata_2ms = ifdata(1:cor_1ms*2);% 截断到2ms
% 求循环相关
for i = 1:Freq_Range/f_delta
    f_temp = f_IF - Freq_Range/2 + i*f_delta + freq_acq;% e5a频率估计值
    DDCData = ifdata_2ms.*exp(1i*2*pi*f_temp/SystemClockFre*(1:cor_1ms*2)).';%  下变频到中频
    cor_sum2 = zeros(Code_Range/code_delta, 1);
    for j = 1:Code_Range/code_delta
        index = ceil((((1:cor_1ms*2)+j*code_delta-Code_Range/2 - code_acq + 1)/cor_1ms*GALILEO_CACODE_NUM));% 将prn码从20460个转换到58000个的转换序列
        index = mod(index -1, GALILEO_CACODE_NUM*2)+1;
        cor_sum2(j) = sum(DDCData.*code_2ms(index));
    end
    cor_sum(i, :) = cor_sum2.';
end
% 求捕获的频率和码相位
cor_pow = cor_sum.*conj(cor_sum);
[r, w] = find(cor_pow == max(cor_pow, [], 'all'));
freq_acqui = - Freq_Range/2 + r*f_delta;% 捕获的多普勒频率
code_phase = w/cor_1ms*GALILEO_CACODE_NUM;% 码相位
end