clc
clear
close all

%% x1 x2 already two imf  
TR = 1;
sliding_window_size_seconds = [30,60,100];
f1 = 1/50; f2 = 1/25;f3 = 1/20;
tf_sec = 300;
[h1] = function_simulation_Zhuang_2019_NI(f1,f2,f3,TR,tf_sec,sliding_window_size_seconds);