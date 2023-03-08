function [] = initial_paper_switch(fft_seg, config, GC_empirical, swtich_vec, t_min_switch, num_switch)
%% add switch setting 
%% use all tapers observatoins rather than using only 1 taper in the generating observation 

fr = config.f_l;
sr = config.sample_r;
t_min_cofig = config.t_min;
t_min_end = t_min_cofig(end);

stp = t_min_end/(length(GC_empirical)-1);
t_min_GCoh = 0: stp: t_min_end;

%% find index of switch
ind_switch_GCoh = zeros(length(t_min_switch), 1);
for i=1: length(t_min_switch)
    abs_diff_i = abs(t_min_GCoh-t_min_switch(i));
    ind_switch_GCoh(i) = find(abs_diff_i == min(abs_diff_i));

end

%% generate Y_k and estimate mu
m_ = size(fft_seg);
num_tpr = m_(2);

m = size(fft_seg{1});
num_ch = m(2);

mm = size(fft_seg{1}{1});
num_chunck = mm(1);
num_win_in_chunk = mm(2);

%% extract data of all tapers 
MU_ch_all_tpr = zeros(num_ch,1, num_tpr);
FFT_ch_all_tpr = zeros(num_ch, num_chunck*num_win_in_chunk, num_tpr);
for tpr_ = 1: num_tpr
    fft_seg_tpr = fft_seg{tpr_};

    fft_ch = zeros(num_ch, num_chunck*num_win_in_chunk);
    mu_ch = zeros(num_ch,1);


    for i=1: num_ch
        temp_fft_ch = fft_seg_tpr{i};

        temp_ch = zeros(1, num_chunck*num_win_in_chunk);

        for j=1: num_chunck
            ind_ = (j-1)*num_win_in_chunk+1: j*num_win_in_chunk;
            temp_ch(1, ind_) = temp_fft_ch(j,:);

        end
        mu_ch(i,1) = mean(temp_ch);
        fft_ch(i, :) = temp_ch;
    end
    MU_ch_all_tpr(:, :, tpr_) = mu_ch;
    FFT_ch_all_tpr(:, :, tpr_) = fft_ch;
end

MU_ch_ave = mean(MU_ch_all_tpr, 3);
FFT_ch_ave = mean(FFT_ch_all_tpr, 3);

%% this part is for switch - calculate switch points 

empirical_switch_pionts = zeros(1, num_chunck*num_win_in_chunk);

for j=1: num_chunck

    if num_switch==3
        if ind_switch_GCoh(1) <=j && j< ind_switch_GCoh(2)
            switch_val = swtich_vec(1);
        elseif ind_switch_GCoh(2) <=j && j< ind_switch_GCoh(3)
            switch_val = swtich_vec(2);
        elseif ind_switch_GCoh(3) <=j && j< ind_switch_GCoh(4)
            switch_val = swtich_vec(3);
        elseif ind_switch_GCoh(4) <=j && j< ind_switch_GCoh(5)
            switch_val = swtich_vec(4);
        elseif ind_switch_GCoh(5) <=j && j< ind_switch_GCoh(6)
            switch_val = swtich_vec(5);
        end
    elseif num_switch ==2
        if ind_switch_GCoh(1) <=j && j< ind_switch_GCoh(2)
            switch_val = swtich_vec(1);
        elseif ind_switch_GCoh(2) <=j && j< ind_switch_GCoh(3)
            switch_val = swtich_vec(2);
        elseif ind_switch_GCoh(3) <=j && j< ind_switch_GCoh(4)
            switch_val = swtich_vec(3);
        end

    elseif num_switch==1 % have 1 switch - SS-Coh

        if ind_switch_GCoh(1) <=j && j< ind_switch_GCoh(2)
            switch_val = swtich_vec(1);
        end
    end

    empirical_switch_pionts(1,(j-1)*num_win_in_chunk+1:j*num_win_in_chunk) = switch_val;
end


figure()
plot(empirical_switch_pionts)
ylim([0, 4])
xlim([0, length(empirical_switch_pionts)])
xlabel('sample')
ylabel('Switch')

%% end of switch calculation of switch points 


Y_k = FFT_ch_ave.';

mu_init = MU_ch_ave;

%% here I can find index of the input for each switch
fft_ch_switch = cell(1, num_switch);
Y_k_switch = cell(1, num_switch);
D_init_switch = cell(1, num_switch);
L_init_switch = cell(1, num_switch);

for k=1: num_switch
    ind_k = find(empirical_switch_pionts==k);
    fft_ch_k = FFT_ch_ave(:, ind_k);
    fft_ch_switch{1, k} = fft_ch_k;
    Y_k_switch{1, k} = fft_ch_k.';
end


%
% %%  cross spectral

for k=1: num_switch
    fft_ch_k = fft_ch_switch{1, k};

    cross_spect_mat_ch = zeros(num_tpr, num_tpr);
    for i=1:num_tpr
        for j=1:num_tpr

            X_i = fft_ch_k(i,:);
            X_j = fft_ch_k(j,:);
            K = length(X_i);

            cross_spect_mat_ch(i,j) = (1/K)*sum(X_i.*conj(X_j));
        end
    end

    %% L, D
    [L_ch , D_ch] = eig(cross_spect_mat_ch);


    % sort
    D_sort_ch = zeros(num_tpr);
    L_sort_ch = zeros(num_tpr);
    for i=1 : num_tpr
        D_sort_ch(i , i) = D_ch(end - (i-1) , end - (i-1));
        L_sort_ch(: , i) = L_ch(: , end - (i-1));
    end
    D_init_ch = D_sort_ch;
    L_init_ch = L_sort_ch;

    L_init_switch{1, k} = L_init_ch;
    D_init_switch{1, k} = D_init_ch;

end



cfg_Init.Y_k = Y_k;
cfg_Init.GC_empirical = GC_empirical;


cfg_Init.mu = mu_init;
cfg_Init.L = L_init_switch;
cfg_Init.D = D_init_switch;
% cfg_Init.incongruent_vec = incongruent_vec;
% cfg_Init.a = a_init;
% cfg_Init.b = b_init;
% cfg_Init.x0 = x0;
cfg_Init.empirical_switch_pionts = empirical_switch_pionts;
cfg_Init.cfgGen = config;



str_save_yk = sprintf('Anethesia_switch_%d_Yk_chNum%d_fr%d_sr%d_winSec%d_segNum%d_TRIAL.mat',num_switch, config.ch_num, fr, sr, config.win_sec, config.seg_num);% using 2 tasks
save(str_save_yk, 'Y_k')

str_save_cofig = sprintf('Anethesia_switch_%d_cfg_Init_chNum%d_fr%d_sr%d_winSec%d_segNum%d_TRIAL.mat',num_switch, config.ch_num, fr, sr, config.win_sec, config.seg_num); % using 2 tasks
save(str_save_cofig, 'cfg_Init');

end