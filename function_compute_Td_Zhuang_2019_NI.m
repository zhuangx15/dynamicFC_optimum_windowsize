function [window_size] = function_compute_Td_Zhuang_2019_NI(tc,TR,max_imf)
%%=========================================================================
% tc: time series;
% window_size: in TR
% IMPORTANT: function emd is not included and is from Patrick Flandrin, 
%            which can be downloaded at http://perso.ens-lyon.fr/patrick.flandrin/emd.html
%%============================================================================

    [tdim,Nq] = size(tc);
    tc = zscore(tc);
    td_t_in_s = zeros(Nq,tdim);
    for j = 1:Nq
        %     j
        tc_1_tmp = function_energy_normalize(tc(:,j));
        tc_1_tmp_flip = flip(tc_1_tmp);
        tc_1 = [tc_1_tmp_flip;tc_1_tmp;tc_1_tmp_flip];
        tdim_f = size(tc_1,1);
        if sum(isnan(tc_1))>0
            continue;
        end
        [energy_tc1,~,~,~,imf_1] = function_calculate_mean_period_mean_energy_V2(tc_1,max_imf,TR,'hilbert');
        N_imf = size(imf_1,1);
        energy_tc1 = energy_tc1./sum(energy_tc1);
        pe1 = zeros(N_imf,tdim_f);
        for ii = 1:N_imf
            if ii<=size(imf_1,1)
                [~,pe1(ii,1:tdim_f-1)] = function_calculate_instantaneous_frequency(imf_1(ii,:),TR);
            end
        end
        td_in_s = energy_tc1'*pe1;
        td_t_in_s(j,:) = td_in_s(tdim+1:tdim*2);
    end
    window_size = cell(Nq,Nq);
    for j = 1:Nq
        td_in_s1 = td_t_in_s(j,:);
        for k = 1:j-1
            td_in_s2 = td_t_in_s(k,:);
            td_in_s_tmp = [td_in_s1;td_in_s2];
            td = (max(td_in_s_tmp)./TR);  %in TR
            window_size{j,k} = td;  % in TR
        end
    end
end

function [mean_energy,mean_period,log_mean_energy,log_mean_period,imf_TC] = ...
    function_calculate_mean_period_mean_energy_V2(TCEN,N_imf_each_ica,TR,type)
    if nargin < 4
        type = 'hilbert';
    end
    TCEN = TCEN(:);
    imf_TC = emd(TCEN);
    N1 = size(imf_TC,1);
    tdim_t = size(TCEN,1);
    mean_energy = zeros(N1,1);
    mean_period = zeros(N1,1);
    for j = 1:N1
        if j>N1
            break;
        end
        imf_tmp = imf_TC(j,:);
        mean_energy(j,1) = sum(imf_tmp.^2)/tdim_t;
        if strcmp(type,'peaks') == 1
            num_peaks = length(findpeaks(imf_tmp));
            mean_period(j,1) = tdim_t/num_peaks*TR;
        elseif strcmp(type,'fft') == 1
            [P_freq,freq] = fourier_transform_plot_f_vs_amplitude(imf_tmp,TR,tdim_t,0);
            P_freq_sum = trapz(freq,P_freq);
            P_freq = P_freq./P_freq_sum;
            period = 1./freq;
            P_period_times_period = - P_freq./ (period);
            mean_period(j,1) = trapz(period(2:end),P_period_times_period(2:end));
        elseif strcmp(type,'hilbert') == 1
            P_hilbert = hilbert(imf_tmp);
            freq = abs(1/TR/(2*pi)*diff(unwrap(angle(P_hilbert)))); %instanteous frequency
            %         freq = sort(freq);
            period = 1./freq;
            mean_period(j,1) = mean(period);
        end
    end
    if N1<N_imf_each_ica
        N_imf_each_ica_f = N1;
    else
        cum_mean_energy = cumsum(mean_energy);
        per_cum_mean_energy = cum_mean_energy ./ cum_mean_energy(end);
        ind = find(per_cum_mean_energy>0.99);
        N_imf_each_ica_f = max([N_imf_each_ica,ind(1)]);
    end
    mean_energy = mean_energy(1:N_imf_each_ica_f,:);
    mean_period = mean_period(1:N_imf_each_ica_f,:);
    imf_TC = imf_TC(1:N_imf_each_ica_f,:);
    log_mean_energy = log(mean_energy);
    log_mean_period = log(mean_period);
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

