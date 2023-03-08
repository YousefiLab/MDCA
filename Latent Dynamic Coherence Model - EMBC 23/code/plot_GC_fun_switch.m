function plot_GC_fun_switch(GC , config, fnt_size, t_min_switch)
% PLOT_GC_FUN(GC , config , method_GC) function plot GC 
% if we have an interval for desired freq our plot is 2D
% if we only want to plot an individual freq, we have a 1D plot
%INPUTS
%  GC a cell that contains GC for different frequencies and window times
%
% config is a structure with 9 fields in it
% 1 - ch_num :    number of channels
% 2 - Fs :        sampling rate of data
% 3 - sample_r :  chosen sample rate for making windows and segments (this is not sample rate of data)
% 4 - win_sec :   this value multiply to sample_r = length of window uses for chunking data
% 5 - win_length  (sample) length of window that we want to calc Global Coherence on them
% 6 - seg_num     (number) number of segments that we want to devide a window into them  (number)
% 7 - seg_length  (sample) number of samples of segments
% 8 - f_l         lower bound of desired freq interval
% 9 - f_u         upper bound of desired freq interval
% 10 - method_GC  specify the method that we use
% with method_GC, the method that we want use to calculate GC is selected
% there are two options,
% first = Proposed, means without subtracting from mean FFT
% second = PNAS, based PNAS Brown paper and subtracting FFT results from FFT's mean
% 11 - over_lap : (of 1) length of overlap window based win length 


%%%% Extract config info
ch_num = config.ch_num;
sample_r = config.sample_r;
win_length = config.win_length;
seg_num = config.seg_num;
f_l = config.f_l;
f_u = config.f_u;
over_lap = config.over_lap;

m = size(GC);
%%% m(1) = num of desired frequencies
%%% m(2) = num of windows

f_ind = [f_l f_u];

win_sec = win_length./sample_r;
% data_length = m(2)*win_sec;
% t_ind = win_sec/2 :win_sec: data_length-win_sec/2;
% % t_ind = t_ind./(60*60); %%% change second scale to hour
% t_ind = t_ind./(60); %%% change second scale to min


t_min = config.t_min;
step = t_min(end)/length(GC);
t_ind = 0 :step : t_min(end)-step;


%%% 2d plot
if m(1)>1
    
    figure
%     figure('units','normalized','outerposition',[0 0 1 1]),


    

    imagesc(f_ind,t_ind,GC.'),colorbar;
    str_tit_1 = sprintf('Global Coherence for %d channels',ch_num);
    str_tit_2 = sprintf('window length = %d sec * number of segment per each window = %d * overlap = %.0f%%',...
        win_sec , seg_num , over_lap*100);
    xlabel('Freq (Hz)'),
    ylabel('Time (mins)')
    ax = gca;
    % xtick
    x_num_step = 3;
    ylim([0 t_ind(end)])
    y_tick_values = fun_finding_tick(ax, x_num_step, 'y'); % nested fun
%     ax.YTick = y_tick_values;
    
%     ax.XTick = [10 20];
    set(gca,'FontSize',fnt_size)
    
    
%     title({str_tit_1 , str_tit_2})
    view(-90,90);
    colormap('jet')
    %% % 1d plot
else
    figure
%     figure('units','normalized','outerposition',[0 0 1 1]);
    plot(t_ind , GC , 'LineWidth',4 , 'color' , 'b')
    xlim([0 t_ind(end)]) , ylim([0 1])
    xlabel('Time (mins)'),ylabel('GC');
    
    % plot vertical lines 
    for ii=1: length(t_min_switch)
        hold on
        xline(t_min_switch(ii))
    end

     set(gca,'FontSize',fnt_size)
    ax = gca;
    % xtick
    x_num_step = 3;
    x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
    ax.XTick = x_tick_values;
    
    str_tit_1 = sprintf('GC');
    str_tit_2 = sprintf('%d channels at freq = %d * window length = %d sec * segments per each window = %d * overlap = %.0f%%'...
        ,ch_num, f_l ,win_sec , seg_num , over_lap*100);
%     title({str_tit_1 , str_tit_2});

%     legend('Original Signal','Smoothed Signal')
    
    
end

set(gcf, 'PaperPosition', [0 0 6 4]);

print('GC_epirical','-dpng','-r600')

end



%% nested functions 
% finding the proper val for ticks
function tick_val = fun_finding_tick(ax, num_step, axis)

switch axis
    case 'x'
        lim_val = ax.XLim;
    case 'y'
        lim_val = ax.YLim;
end

step_tick = ((lim_val(end) - lim_val(1))./num_step);
cnt = find_order(step_tick); % nested fun
if cnt> 1
    step_tick_n = round(step_tick, -(cnt-1));
    tick_val = round(lim_val(1) + (step_tick_n./2), -(cnt-1)) : step_tick_n: lim_val(end);
else
    step_tick_n = round(step_tick, cnt);
    tick_val = round(lim_val(1) + (step_tick_n./2), cnt) : step_tick_n: lim_val(end);
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
