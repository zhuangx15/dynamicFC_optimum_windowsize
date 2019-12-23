function [h1] = function_simulation_Zhuang_2019_NI(f1,f2,f3,TR,tf_sec,fixed_window_size_seconds)
%% generate simulated sinusoid singals with frequency shifts as in Zhuang et al., 2019 NI (detials in section 2.2);

lw = 2;  % lineWidth in plot  
fontSize = 16;
sliding_window_size_vec = ceil(fixed_window_size_seconds/TR);
N_win = numel(sliding_window_size_vec);
t = [1:TR:tf_sec]'; tdim = size(t,1);
x1 = zeros(tdim,1);
index1 = 1:floor(200/TR);
x1(index1) = cos(2*pi*f1*t(index1));
index2 = floor(200/TR)+1:tdim;
x1(index2) = cos(2*pi*f3*t(index2));
x2 = cos(2*pi*f2*t);
[x1] = function_energy_normalize(x1);
[x2] = function_energy_normalize(x2);
h1 = figure(100); set(h1,'position',[50,50,2400,1000]);
figure(h1);subplot(221),plot(t,x1,t,x2,'LineWidth',lw);  grid on;
legend({'y^1','y^2'},'fontSize',12);
title(['Time Series: TR=',num2str(TR),'s'],'FontSize',fontSize);  xlabel('s','FontSize',fontSize);
ylabel('Intensity','fontSize',fontSize);
static_corr = corr(x1,x2);

[Px] = fourier_transform_plot_f_vs_amplitude(x1,TR,tdim);
[Py,f] = fourier_transform_plot_f_vs_amplitude(x2,TR,tdim);
figure(h1),subplot(222),plot(f,Px,f,Py,'LineWidth',lw); grid on;
legend({'y^1','y^2'},'fontSize',12);
title('Frequency Domain','FontSize',fontSize);  xlabel('Hz','FontSize',fontSize);
ylabel('Intensity','fontSize',fontSize);

x11 = [flip(x1);x1;flip(x1)];
x22 = [flip(x2);x2;flip(x2)];
[x11] = function_energy_normalize(x11);
[x22] = function_energy_normalize(x22);
fx = function_calculate_instantaneous_frequency(x11,1);
pe1 = 1./(fx(tdim+1:tdim*2));
fy = function_calculate_instantaneous_frequency(x22,1);
pe2 = 1./fy(tdim+1:tdim*2);
td= round(max([pe1';pe2'])); td = td(:); %in TR;
td(td>(tf_sec/TR)) = floor(tf_sec/TR);
figure(h1);subplot(224),plot(t,pe1,t,pe2,'LineWidth',lw); hold on; plot(t,td,'g--','LineWidth',lw+1);
grid on; %hold off;
xlabel('s','fontSize',fontSize);
ylabel('TR','fontsize',fontSize);
legend_display = cell(N_win,1);
for j = 1:N_win
    plot(t(1:tdim-1),sliding_window_size_vec(j)*ones(tdim-1,1),'--','LineWidth',lw);
    legend_display{j,1} = ['Window-Size:',num2str(fixed_window_size_seconds(j)),'s'];
end
hold off; ylabel('TR','FontSize',fontSize); title('Window-Size','FontSize',fontSize);
legend(['instantaneous period y^1';'instantaneous period y^2';'Window-Size: optimum';legend_display],'fontSize',12);
%% correlation calculate using instantaneous perios as the optimal window size
opt_window_corr = NaN(tdim,1);
for i = 1:tdim-1
    if round(i-td(i,1)/2) >= 1 && round(i+td(i,1)/2) <= tdim
        opt_window_corr(i,1) = corr(x1(round(i-td(i,1)/2):round(i+td(i,1)/2)),x2(round(i-td(i,1)/2):round(i+td(i,1)/2)));
    end
end
figure(h1);subplot(223);
h(1) = plot(t,nan(tdim,1)); hold on;
h(2) = plot(t,nan(tdim,1)); hold on; h(3) = plot(t,opt_window_corr,'g','LineWidth',lw); hold on;
%% sliding window
for j = 1:N_win
    sliding_window_corr = NaN(tdim,1);
    wn = sliding_window_size_vec(j);
    for i = 1:tdim
        if round(i-wn/2) >= 1 && round(i+wn/2) <= tdim
            sliding_window_corr(i,j) = corr(x1(round(i-wn/2):round(i+wn/2)),x2(round(i-wn/2):round(i+wn/2)));
        end
    end
    h(3+j) = plot(t,sliding_window_corr(:,j),'--','LineWidth',lw);
end
h(N_win+3+1)= plot(t,static_corr*ones(tdim,1),'m--','LineWidth',1.25);
hold off;
legend(h(3:end),['Window Size: Optimum';legend_display;'Static Correlation'],'fontSize',12);
grid on;
title('\fontsize{16}Dynamic FC value');
xlabel('s','fontSize',fontSize);
ylabel('Correlation','fontSize',fontSize);
set(gcf,'color','w');
end

function [TCEN,energy] = function_energy_normalize(TC)
% tc should be a vector
% mean(TCEN.^2) = 1;
TC = TC(:);
tdim = size(TC,1);
TCEN = TC./sqrt(sum(TC.^2)./tdim);
energy = 1/tdim*sum(TCEN.^2);
end

function [f,period] = function_calculate_instantaneous_frequency(x1,TR)
% energy normalized x1;
P1_hilbert = hilbert(x1);
tdim = numel(x1);
f = abs(1/TR/(2*pi)*diff(unwrap(angle(P1_hilbert))));  %in Hz
period = 1./f; % in sec
period(period>(TR*tdim)) = 60; % if the period is greater then the entire TR*tdim;
end


function [P1,f] = fourier_transform_plot_f_vs_amplitude(ft,deltat,L)
ft = ft(:);
if mean(ft) ~= 0
    ft = ft - mean(ft);
end
Fs = 1/deltat;             % sampling frequency
t = (0:L-1)'*deltat;        % Time vector
Fw = fft(ft,L);
P2 = abs(Fw/L);
P1 = P2(1:round(L/2));
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*((0:(ceil(L/2)-1))')/L;
end
