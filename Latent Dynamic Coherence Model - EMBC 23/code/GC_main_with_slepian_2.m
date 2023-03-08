clc
clear
close all

load('E:\oneDrive_WPI\OneDrive - Worcester Polytechnic Institute (wpi.edu)\-0\switching_data\ali_3.mat');
% 
data_det = strct_switch_data.data_det; 
data_det = data_det.';
t_sec = strct_switch_data.t_sec;
Fs = strct_switch_data.Fs;
Ts = strct_switch_data.Ts;

% data_reza = data_reza.data;
data_ = data_det.';
l_d = size(data_det, 2);

t_s = t_sec;
t_min = t_s/60;
ch_n = size(data_det, 2);

% data_det =data_;
method_GC = 'PNAS';

% this method affetcs only on angle_eig_vec_fun
% method_angle = 1;
method_angle = 2;



%% config 
%%% analyze configuration
config.ch_num = ch_n;                                       % number of channels
config.Fs = Fs;                                         % sampling rate of data
config.sample_r = Fs;                                   % size of slepian window
config.win_sec = ch_n*1;
config.win_length = config.win_sec*config.sample_r;      % (sample) length of window that we want to calc Global Coherence on them
config.seg_num = config.win_sec;                                      % (number) number of segments that we want to devide a window into them  (number)
config.seg_length = config.win_length./config.seg_num;   % (sample)
config.f_l = 2;                                          % lower bound of desired freq interval
config.f_u = config.f_l + 28;                            % upper bound of desired freq interval
config.over_lap = .8;                                   % (of 1) length of overlap window based win length
config.method_angle = 2;
config.num_slepian = 1;
config.halfbandwidth = 3.5;
config.t_min = t_min;

%%% data extraction
m=size(data_det);

%%% Time Plot
plt_time = 0;
while plt_time == 1
    plt_time = 0;
    ch_plt_vec = [6 40];
    fnt_size = 25;
    
    for i=1 : config.ch_num
%         ch_plt = ch_plt_vec(i);
        figure('units','normalized','outerposition',[0 0 1 1]),
        plot(t_min, data_det(:,i), 'Color', [27, 27, 28]./255)
        
        ylim([-3000 3000])
        
        set(gca,'FontSize', fnt_size)
        
% 
        set(gca, 'YTIck', [-500:500:500])
        set(gca, 'YAxisLocation', 'right')
        
        


        xlabel('Time (mins)')
        ylabel('Amplitude')
        xlim([0 t_min(end)])
        ylim([-800 800])
        
        set(gca, 'XTIck', [20:40:140])
        
        
        set(gcf, 'PaperPosition', [0 0 6 4]);
%         print('time_40','-dpng','-r600')
        
        
        
    end
    close all
end


data_chs = data_det; 

% plt_time_data(data_chs, t_min)
% config.num_ch = 19;
% normalized time data
% max_val = max(abs(data_chs(:)));
% data_chs = data_chs./max_val;

 %% multitaper & ordinary spectrogram
plt_mltp = 0;
while plt_mltp==1
    plt_mtlp =0;
    % input structure of multitaper spectrogram
    cfg_mtp_spect = struct;
    
    cfg_mtp_spect.T = T;               % time indices
    cfg_mtp_spect.ch_des_num = length(ch_rand);     % desired channel for plotting spectrogram
    cfg_mtp_spect.f_des_l = 1;         % (freq(Hz)) the lower range for plotting spectrogram
    cfg_mtp_spect.f_des_u = 29;         % (freq(Hz)) the lower range for plotting spectrogram
    cfg_mtp_spect.mtp_NW = 3.5;         % NW (halfbandwidth) value for calculating multitaper
    cfg_mtp_spect.mtp_win_length = 60;  % (sec) win length for calculating multitaper
    cfg_mtp_spect.mtp_over_lap = .3;   % (of 1) length of overlap window based win length
    
    for i =1: length(ch_rand)
        cfg_mtp_spect.ch_des_num = ch_rand(i);       % desired channel for plotting spectrogram
        
        % this fuction plot spectrogram based multitaper approach for desired channel
        multiTaper_spectrogram(data_det , cfg_mtp_spect)
    end
    % this fuction plot spectrogram based multitaper approach for desired channel
    
end



%% calc GC and Igen info
% [GC , sorted_eig_info] = calc_GC_steps(data_chs , config);
[GC_overal , sorted_eig_info, fft_seg_mean_removed] = calc_GC_steps_rsf(data_chs , config);

% str_save = sprintf('anesthesia_ch%d_win%d_overLap%.1f_fL%d_fU%d_numSlep%d.mat', config.ch_num, config.win_sec, config.over_lap, config.f_l, config.f_u, config.num_slepian);
% save(str_save, 'sorted_eig_info')

% we can load instead of calculating again GC adn EIGEN info if we saved it
% for desired freq before

% struct_sub6_gc = [];

% load('anesthesia_ch32_win8_overLap0_fL11_fU11_numSlep5')
% load('ch32_win32_seg32_frL25_frU25_overLap0')

%% (ANALYZING RESULTS)

win_sec_GCoh_empircal = config.win_sec;
win_sec_observation = config.win_sec / config.seg_num;
fprintf('\n Number of channel: %d' , ch_n);
fprintf(' \n Length of Window for calc GCOh empirical: %d sec', win_sec_GCoh_empircal)
fprintf(' \n Length of subWindow for calc Cov mat of GCOh empirical (Window for SS-Coh): %d sec', win_sec_observation)
fprintf('\n------------------------------\n')

fnt_size = 20;

%%% Global Coherence PLOTTING
plot_GC_fun(GC_overal , config, fnt_size);
% close all;

save('E:\oneDrive_WPI\OneDrive - Worcester Polytechnic Institute (wpi.edu)\-0\figure_1\gc_empr.mat', 'GC_overal', 'config');

%%% calculating angle between largest corresponding eigenvectors over time (1d)
angle_eig_vec_fun(sorted_eig_info , config, fnt_size);
% close all;



%% reza for whole cross spectral/Yk-observation- (PAPER EMBC)
initial_paper(fft_seg_mean_removed{1}, config, GC_overal)
k=1
%% sum eigenvalues
% [max_eig_val , sum_eig_val] = max_sum_eig_value(sorted_eig_info , config );
% close all

% ang_diff_phase_ev(sorted_eig_info , config, fnt_size)
% close all 

% trace_ev(sorted_eig_info , config, fnt_size )
% close all


%% nested functions
function plt_time_data(data, t_min)

my = size(data);
x = t_min;

max_val = max(data(:));
y_lim = [-1 max_val+1];

tick_num_step = 3;
% y_ticks = [0 4];
%% setup values
left_pos = .15; 
bot_pos = .15;
width_val = .7;
height_val = .20;


bot_diff = .22;

line_width = .2;
fnt_size= 14;

%% plots 
h =figure;
ax = subplot(3,1,3);
plot(x, data(:, 1), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 0*bot_diff) width_val height_val];
ylabel('Ch 1')
xlabel('Time (min)')
ylim([-600 600])
xlim([0 t_min(end)])

ax.FontSize = fnt_size;
% xtick
x_num_step = tick_num_step;
x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
ax.XTick = x_tick_values;

% ytick
% y_num_step = tick_num_step-1;
% y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
% ax.YTick = y_tick_values;
ax.YTick = [-300 300];

ax.XAxis.FontSize = fnt_size + 3;
ax.XLabel.FontSize = fnt_size + 3;


ax = subplot(3,1,2);
plot(x, data(:, 4), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 1*bot_diff) width_val height_val];
ylabel('Ch 4')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
% y_num_step = tick_num_step-1;
% y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
% ax.YTick = y_tick_values;

ylim([-600 600])
xlim([0 t_min(end)])
ax.FontSize = fnt_size;
ax.YTick = [-300 300];

ax = subplot(3,1,1);
plot(x, data(:, 3), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 2*bot_diff) width_val height_val];
ylabel('Ch 3')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
y_num_step = tick_num_step-1;
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
ax.YTick = y_tick_values;

% ylim(y_lim)
xlim([0 t_min(end)])
ax.FontSize = fnt_size;

ax = subplot(5,1,2);
plot(x, data(:, 4), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 3*bot_diff) width_val height_val];
ylabel('Ch 4')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
y_num_step = tick_num_step-1;
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
ax.YTick = y_tick_values;

% ylim(y_lim)
xlim([0 t_min(end)])
ax.FontSize = fnt_size;

ax = subplot(5,1,1);
plot(x, data(:, 5), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 4*bot_diff) width_val height_val];
ylabel('Ch 5')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
y_num_step = tick_num_step-1;
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
ax.YTick = y_tick_values;

% ylim(y_lim)
xlim([0 t_min(end)])
ax.FontSize = fnt_size;



% [left bottom width height]

end

% finding the proper val for ticks
function tick_val = fun_finding_tick(ax, num_step, axis)

switch axis
    case 'x'
        lim_val = ax.XLim;
        div_num = 2;
    case 'y'
        lim_val = ax.YLim;
        div_num = 3;
end

step_tick = ((lim_val(end) - lim_val(1))./num_step);
cnt = find_order(step_tick); % nested fun
if cnt> 1
    step_tick_n = round(step_tick, -(cnt-1));
    tick_val = round(lim_val(1) + (step_tick_n./div_num), -(cnt-1)) : step_tick_n: lim_val(end);
else
    step_tick_n = round(step_tick, cnt);
    tick_val = round(lim_val(1) + (step_tick_n./div_num), cnt) : step_tick_n: lim_val(end);
end


end

% finding the order of input value
function cnt = find_order(x)

if x> 1
    stp_con = 0;
    cnt = -1;
    while stp_con == 0
        cnt = cnt +1;
        a = floor(x./10^cnt);
        
        if a == 0
            stp_con=1;
        end
        
    end
else
    stp_con = 0;
    cnt = 0;
    while stp_con == 0
        cnt = cnt +1;
        a = floor(x./10^(-cnt));
        
        if a > 0
            stp_con=1;
        end
        
    end
end

end
