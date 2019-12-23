clc
clear

load test_tc.mat tc;
TR = 0.72;
max_imf = 10;
tdim = size(tc,1);
[window_size] = function_compute_Td_Zhuang_2019_NI(tc,TR,max_imf);