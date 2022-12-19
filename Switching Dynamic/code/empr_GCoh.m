clc
clear
close all



% mu1 = [1 2];
% Sigma1 = [2 0; 0 0.5];
% mu2 = [-3 -5];
% Sigma2 = [1 0;0 1];
% rng(1); % For reproducibility
% X = [mvnrnd(mu1,Sigma1,1000);mvnrnd(mu2,Sigma2,1000)];




% str_data = '8_12_2min_11973';  % not good 
% str_data = 'ali_reading_6_16_8529'; % not good
% str_data = 'ali_reading_6_16_eye_8526'; % good - 16 channel
% str_data = 'reza_3';
% str_data = 'ali_0001';
% str_data = 'ali_11_30_22';
% str_data = 'reza2_12_2_22';
str_data = 'ken';


% str_data = '8_12_2min_4113';
% str_data = '8_12_2min_4699';
% str_data = '8_12_2min_4707';
% str_data = 'A1';
% str_data = 'A2';

load(str_data)

data_det = strct_switch_data.data_det;
t_sec = strct_switch_data.t_sec;
Fs = strct_switch_data.Fs;
Ts = strct_switch_data.Ts;

t_min = t_sec/60;


ch_vec = 1: size(data_det, 1);
ch_vec2 = 1: 2*size(data_det, 1);
%% 
disp('Calculate GCoh emprical')
data_det = data_det.';
win = 8;
params = [];
bw = 5;
tpr_num = 7;
params.tapers = [bw tpr_num];
params.Fs = Fs;
params.fpass = [2 30];

win_length = win;
win_overlap = 0.8;
ind_win = calc_sec_intrvl(t_sec, win_length, win_overlap);

data_ = data_det(ind_win(1, 1): ind_win(2, 1), :);
[Sc, Cmat, Ctot, Cvec, Cent, f]=CrossSpecMatc(data_ ,win, params);

f_win = f;
l_win = length(ind_win);
if win_overlap ~= 0.0
    l_gcoh = l_win - 5;
else
    l_gcoh = l_win;
end 
l_freq = length(f_win);
t_win = zeros(1, l_gcoh);
gcoh_win = zeros(l_gcoh, l_freq);

num_ch = size(data_det, 2);
eVec_fr_l = zeros(l_gcoh, num_ch); 
eVec_fr_u = zeros(l_gcoh, num_ch); 
eVec_fr_l_u = zeros(l_gcoh, 2*num_ch); 
fr_l = 9; 
fr_u = 25; 

diff_fr_l = abs(f_win - fr_l);
ind_fr_l = find(diff_fr_l == min(diff_fr_l));
f_l = f_win(ind_fr_l);
diff_fr_u = abs(f_win - fr_u);
ind_fr_u = find(diff_fr_u == min(diff_fr_u));
f_u = f_win(ind_fr_u);

for i=1: l_gcoh
    str_progress = sprintf('%d/%d', i, l_win);
    disp(str_progress);
    
    ind_l = ind_win(1, i);
    ind_u = ind_win(2, i);
    
    t_l = t_sec(ind_l);
    t_u = t_sec(ind_u);
    t_avg = mean([t_l t_u]);
    t_win(1, i) = t_avg; 
    
    data_win = data_det(ind_l: ind_u-1, :);
    
%     for j=1: 16
%         
%        figure 
%        plot(data_win(:, j))++
%     end
    
    % data (in form samples x channels x trials) 
    % win (duration of non-overlapping window)
    % params: structure with fields tapers, pad, Fs, fpass
    [Sc, Cmat, Ctot, Cvec, Cent, f]=CrossSpecMatc(data_win, win, params);
    % Sc = (cross spectral matrix frequency x channels x channels)
    % Cmat Coherence matrix frequency x channels x channels
    % Ctot Total coherence SV(1)^2/sum(SV^2) (frequency)
    % Cvec leading Eigenvector (frequency x channels)
    % Cent A different measure of total coherence: GM/AM of SV^2s f (frequencies)
    gcoh_win(i, :) = Ctot;
    
    
    eVec_fr_l_i = abs(Cvec(ind_fr_l, :));
    eVec_fr_u_i = abs(Cvec(ind_fr_u, :));
    
    eVec_fr_l(i, :) = eVec_fr_l_i;
    eVec_fr_u(i, :) = eVec_fr_u_i;
    
    eVec_fr_l_u_i = [eVec_fr_l_i, eVec_fr_u_i];
    eVec_fr_l_u(i, :) = eVec_fr_l_u_i;
    
        
end


plt_empr_Gcoh_freq = 1;
if plt_empr_Gcoh_freq
    gcoh_l = gcoh_win(:, ind_fr_l);
    gcoh_u = gcoh_win(:, ind_fr_u);
    
    figure()
    plot(t_win/60, gcoh_l, 'LineWidth', 2);
    xlabel('Time')
    str_tit = sprintf('Gcoh - %.1f Hz', f_l);
    title(str_tit)
    
    figure()
    plot(t_win/60, gcoh_u, 'LineWidth', 2);
    xlabel('Time')
    str_tit = sprintf('Gcoh - %.1f Hz', f_u);
    title(str_tit)
end

% str_eVec = [];
% 
% str_eVec.eVec_10 = eVec_fr_l;
% str_eVec.eVec_25 = eVec_fr_u;
% str_eVec.fr_l = f_l;
% str_eVec.fr_u = f_u;
% str_eVec.t_min = t_win/60;
% 
% save('str_eVec_swtichData.mat', 'str_eVec')


plt_ = 1;
if plt_==1
    %% plot GCoh
    figure()
    imagesc(t_win/60, f_win, (gcoh_win.'))
    xlabel('Time (min)')
    ylabel('Freq (Hz)')
    colormap jet
    colorbar()
    set(gca,'YDir','normal')
    set(gca,'FontSize',20)
    
    %% plot leading eVec for fr l
    figure()
    imagesc(t_win/60, ch_vec, (eVec_fr_l.'))
    xlabel('Time (min)')
    ylabel('Ch')
    colormap jet
    colorbar()
    set(gca,'YDir','normal')
    str_tit = sprintf('lead eVec (%.2f)', f_l);
    title(str_tit)
    set(gca,'FontSize',20)
    
    %% plot leading eVec for fr u
    figure()
    imagesc(t_win/60, ch_vec, (eVec_fr_u.'))
    xlabel('Time (min)')
    ylabel('Ch')
    colormap jet
    colorbar()
    set(gca,'YDir','normal')
    str_tit = sprintf('lead eVec (%.2f)', f_u);
    title(str_tit)
    set(gca,'FontSize',20)
    
    %% plot leading eVec for fr l, u
    figure()
    imagesc(t_win/60, ch_vec2, (eVec_fr_l_u.'))
    xlabel('Time (min)')
    ylabel('Ch')
    colormap jet
    colorbar()
    set(gca,'YDir','normal')
    str_tit = sprintf('lead eVec (%.1f-%.1f)', f_l, f_u);
    title(str_tit)
    set(gca,'FontSize',20)
end
%% save GMM data

dict_dat = [];
dict_dat.eVec_fr_l_u = eVec_fr_l_u;
dict_dat.t_min = t_win/60;
dict_dat.fl_fu = [f_l, f_u];
str_save = sprintf('%s_gmm.mat', str_data);
save(str_save, 'dict_dat')


u=1; 

%% nested functions
function [ind_win] = calc_sec_intrvl(t_sec, win_length, win_overlap)
%% this function calculate index of time based on win length and overlap size of windows 

ind_win = [0;0];

val_strt = 0;
stop_ = 0;
while stop_==0
    diff_strt = abs(t_sec-val_strt);
    ind_l = find(min(diff_strt) == diff_strt);
    
    val_end = val_strt + win_length;
    diff_end = abs(t_sec-val_end);
    ind_u = find(min(diff_end) == diff_end);
    
    if ind_l==ind_u
        stop_ = 1;
    else
        tmp_ = [ind_l; ind_u];
        ind_win = [ind_win, tmp_];
        val_strt = val_strt + (1-win_overlap)*win_length;
    end
    
end
m = size(ind_win);
ind_win = ind_win(:, 2:end-1);
end
