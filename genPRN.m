% @file brief: 生成伪随机码
%% 测试
% bins = hex2bins('acdaee124', 24)
% bin1s = oct2bins('40503', 14)
%%
% array_reduction([1 1 1 1 0 1 1 1 0 1], @xor)
%                               
% %% 生成PRN 读取中频数据
% clear
% [codeai, codeaq, codebi, codebq] = genPRN(7);
% %%
% ifdata = ReadIF('./20210428_A40ms.txt');
% %% 捕获
% SystemClockFre = 5.8e7; %记录回放仪采样率58MHz
% GALILEO_FREQ_COMPENSATE = -0.56439e6; %伽利略的中频频率
% GALILEO_CACODE_FREQ = 10.23e6; %伽利略的CA码速率10.23MHz
% GALILEO_CACODE_NUM = 10230; %伽利略的CA码长度
% 
% f0 = GALILEO_FREQ_COMPENSATE;% 估计的中频频率
% f0_delta = 500;
% fs = 1.5*GALILEO_CACODE_FREQ;% 估计的副载波频率
% fs_delta = 5;
% cor_sum = zeros(1000, 58000);
% 
% % for i = 1:80
% %     f_temp = f0 - fs - 10000 + i*250;% e5a频率估计值
% %     DDCData = ifdata(1:58000).*exp(-1i*2*pi*f_temp/SystemClockFre*(1:58000)).';
% %     parfor j = 1:10230*2
% %         index = ceil(((1:58000)/58000*10230)+0.5*j);
% %         index = mod(index-1, 10230) + 1;
% %         cor_sum(i,j) = sum(codeaq(index).*DDCData);
% %     end
% % end
% index = ceil(((1:58000)/58000*10230));
% parfor i = 1:1000
%     f_temp = f0 - fs - 10000 + i*5;% e5a频率估计值
%     DDCData = ifdata(1:58000).*exp(-1i*2*pi*f_temp/SystemClockFre*(1:58000)).';
%     
%     DDCData_F = fft(DDCData);
%     codeaq_F = fft(codeaq(index));
%     cor_sum(i,:) = (ifft(conj(DDCData_F).*codeaq_F));
% end
% %%
% cor_sum_abs = abs(cor_sum(:));
% %%
% plot(cor_sum.');
% %%
% surf(cor_sum(450:1:500,3.53e4:3.54e4).*conj(cor_sum(450:1:500,3.53e4:3.54e4)),'LineStyle', 'none')
% %%
% surf(cor_sum(1:20:end,1:2:end).*conj(cor_sum(1:20:end,1:2:end)),'LineStyle', 'none')
% %% 
% cor_pow = cor_sum.*conj(cor_sum);
% image(cor_pow.', 'CDataMapping','scaled');
% colorbar
% [r, w] = find(cor_pow == max(cor_pow, [], 'all'));
% %%
% save('./data/cor_sum.mat', 'cor_sum');
%%
function [codeai, codeaq, codebi, codebq] = genPRN(svid)
% @brief: 输出周期为10230的E5信号4路主码
% @param:
% svid - 卫星编号 1~32
% @return:
% codeai - E5a I路PRN码 10230x1
% codeaq - E5a Q路PRN码 10230x1
% codebi - E5b I路PRN码 10230x1
% codebq - E5b Q路PRN码 10230x1
% @warning
StartValueE5a_I ={'30305' '14234' '27213' '20577' '23312' '33463' '15614' '12537' '01527' '30236' ...
    '27344' '07272' '36377' '17046' '06434' '15405' '24252' '11631' '24776' '00630' ...
    '11560' '17272' '27445' '31702' '13012' '14401' '34727' '22627' '30623' '27256' ...
    '01520' '14211' '31465' '22164' '33516' '02737' '21316' '35425' '35633' '24655' ...
    '14054' '27027' '06604' '31455' '34465' '25273' '20763' '31721' '17312' '13277' };
StartValueE5a_Q = {'25652' '05142' '24723' '31751' '27366' '24660' '33655' '27450' '07626' '01705' ...
    '12717' '32122' '16075' '16644' '37556' '02477' '02265' '06430' '25046' '12735' ...
    '04262' '11230' '00037' '06137' '04312' '20606' '11162' '22252' '30533' '24614' ...
    '07767' '32705' '05052' '27553' '03711' '02041' '34775' '05274' '37356' '16205' ...
    '36270' '06600' '26773' '17375' '35267' '36255' '12044' '26442' '21621' '25411' };
StartValueE5b_I = {'07220' '26047' '00252' '17166' '14161' '02540' '01537' '26023' '01725' '20637' ...
    '02364' '27731' '30640' '34174' '06464' '07676' '32231' '10353' '00755' '26077' ...
    '11644' '11537' '35115' '20452' '34645' '25664' '21403' '32253' '02337' '30777' ...
    '27122' '22377' '36175' '33075' '33151' '13134' '07433' '10216' '35466' '02533' ...
    '05351' '30121' '14010' '32576' '30326' '37433' '26022' '35770' '06670' '12017' };
StartValueE5b_Q = {'03331' '06143' '25322' '23371' '00413' '36235' '17750' '04745' '13005' '37140' ...
    '30155' '20237' '03461' '31662' '27146' '05547' '02456' '30013' '00322' '10761' ...
    '26767' '36004' '30713' '07662' '21610' '20134' '11262' '10706' '34143' '11051' ...
    '25460' '17665' '32354' '21230' '20146' '11362' '37246' '16344' '15034' '25471' ...
    '25646' '22157' '04336' '16356' '04075' '02626' '11706' '37011' '27041' '31024' };

%svid = 7;% 7号卫星

c1 = oct2bins('77777');% reg1初始化序列
c2ai = oct2bins(StartValueE5a_I{svid});% reg2 初始化序列
c2aq = oct2bins(StartValueE5a_Q{svid});% reg2 初始化序列
c2bi = oct2bins(StartValueE5b_I{svid});% reg2 初始化序列
c2bq = oct2bins(StartValueE5b_Q{svid});% reg2 初始化序列

% 初始化序列长度为14 后14位为有效数据
c1ai = c1(2:15);
c1aq = c1(2:15);
c1bi = c1(2:15);
c1bq = c1(2:15);

c2ai = c2ai(2:15);
c2aq = c2aq(2:15);
c2bi = c2bi(2:15);
c2bq = c2bq(2:15);

phaseai = 1;
phaseaq = 1;
phasebi = 1;
phasebq= 1;

codeai = zeros(10230,1);
codeaq = zeros(10230,1);
codebi = zeros(10230,1);
codebq = zeros(10230,1);
% 生成C/A码
for i = 1:10230
    [c1ai, c2ai, phaseai, code] = GenPrimaryCodes(c1ai, c2ai, phaseai, 1);
    codeai(i) = 2*code-1;
    [c1aq, c2aq, phaseaq, code] = GenPrimaryCodes(c1aq, c2aq, phaseaq, 2);
    codeaq(i) = 2*code-1;
    [c1bi, c2bi, phasebi, code] = GenPrimaryCodes(c1bi, c2bi, phasebi, 3);
    codebi(i) = 2*code-1;
    [c1bq, c2bq, phasebq, code] = GenPrimaryCodes(c1bq, c2bq, phasebq, 4);
    codebq(i) = 2*code-1;
end
end
% 主码生成
function [c1, c2, phase, code] = GenPrimaryCodes(c1, c2, phase, component)
% @brief: 输出周期为10230的主码
% @param:
% c1 - 移位寄存器1
% c2 - 移位寄存器2
% phase - 码相位 1~10230
% component - 信号类型 1~4 E5aI E5aQ E5bI E5bQ
% @return:
% c1 - 移位寄存器1 更新
% c2 - 移位寄存器2 更新
% phase - 码相位 1~10230 更新
% code - 输出的主码
% @warning
switch component
    case 1
        a1 = oct2bins('40503');% 移位寄存器1的抽头系数
        a2 = oct2bins('50661');% 移位寄存器2的抽头系数
    case 2
        a1 = oct2bins('40503');% 移位寄存器1的抽头系数
        a2 = oct2bins('50661');% 移位寄存器2的抽头系数
    case 3
        a1 = oct2bins('64021');% 移位寄存器1的抽头系数
        a2 = oct2bins('51445');% 移位寄存器2的抽头系数
    case 4
        a1 = oct2bins('64021');% 移位寄存器1的抽头系数
        a2 = oct2bins('43143');% 移位寄存器2的抽头系数
end
a1 = a1(1:14);
a2 = a2(1:14);


code = xor(c1(1),c2(1));
phase = phase + 1;
c1_new = array_reduction(a1.*c1, @xor);
c2_new = array_reduction(a2.*c2, @xor);
c1 = [c1(2:end) c1_new];% 移位
c2 = [c2(2:end) c2_new];% 移位
end

function [bin] = hex2bin(hex)
% @brief: 16进制转换为2进制
% @param: 
% hex - 输入的字符
% @return: 
% bin - 输出的二进制数组
% @warning

switch(hex)
    case {'A', 'a'}
        bin = [1 0 1 0];
    case {'B', 'b'}
        bin = [1 0 1 1];
    case {'C', 'c'}
        bin = [1 1 0 0];
    case {'D', 'd'}
        bin = [1 1 0 1];
    case {'E', 'e'}
        bin = [1 1 1 0];
    case {'F', 'f'}
        bin = [1 1 1 1];
    case '0'
        bin = [0 0 0 0];
    case '1'
        bin = [0 0 0 1];
    case '2'
        bin = [0 0 1 0];
    case '3'
        bin = [0 0 1 1];
    case '4'
        bin = [0 1 0 0];
    case '5'
        bin = [0 1 0 1];
    case '6'
        bin = [0 1 1 0];
    case '7'
        bin = [0 1 1 1];
    case '8'
        bin = [1 0 0 0];
    case '9'
        bin = [1 0 0 1];
end
end

function [bin] = hex2bins(hex, len)
% @brief: 16进制转换为2进制 
% @param: 
% hex - 输入的字符数组
% len - 截短长度
% @return: 
% bin - 输出的二进制数组
% @warning    
bin= arrayfun(@(x) hex2bin(x), hex, 'UniformOutput', 0);
bin = cell2mat(bin);
bin = bin(1:len);
end

function [bin] = oct2bin(oct)
% @brief: 8进制转换为2进制 
% @param: 
% hex - 输入的字符
% @return: 
% bin - 输出的二进制数组
% @warning    
switch(oct)
    case '0'
        bin = [0 0 0];
    case '1'
        bin = [0 0 1];
    case '2'
        bin = [0 1 0];
    case '3'
        bin = [0 1 1];
    case '4'
        bin = [1 0 0];
    case '5'
        bin = [1 0 1];
    case '6'
        bin = [1 1 0];
    case '7'
        bin = [1 1 1];
end
end

function [bin] = oct2bins(oct)
% @brief: 8进制转换为2进制 
% @param: 
% oct - 输入的字符数组
% len - 截断长度
% @return: 
% bin - 输出的二进制数组
% @warning    
bin= arrayfun(@(x) oct2bin(x), oct, 'UniformOutput', 0);
bin = cell2mat(bin);
end

function [result] = array_reduction(array, fun_h)
% @brief: 对数组进行归约
% @param: 
% array - 输入的数组
% fun_h - 具体归约函数句柄 输入参数数量2 输出参数数量1 如xor or and等
% @return: 
% result - 归约结果
% @warning    
result = array(1);
for i = 2:length(array)
    result = fun_h(result, array(i));
end
end