% 读取中频数据
function txt_data = ReadIF(path)
h_file_txt = fopen(path);
temp1 = textscan(h_file_txt, "%d %d");

txt_data = double(temp1{1})+double(temp1{2})*1i;
end