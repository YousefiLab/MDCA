function [max_eig_val , sum_eig_val] = max_sum_eig_value(eig_info , config )
%this function calculate and plot max and summation of all eigenvlues
%
%INPUTS
% 1- eig_info 
% eig_info is a sorted  cell with size 2
% first cell of it includes EIGENVALUES
% second cell of it includes EIGENVECTORS
%
% 2-config
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


%%%% Extract config info
ch_num = config.ch_num;
sample_r = config.sample_r;
win_length = config.win_length;
seg_num = config.seg_num;
f_l = config.f_l;
f_u = config.f_u;
method_GC = config.method_GC; 

eig_val = eig_info{1}; 

m = size(eig_val);
%%% m(1) = number of desired freq num 
%%% m(2) = number of window numbers

win_sec = win_length./sample_r;
% data_length = m(2)*win_sec;
% t_ind = win_sec/2 :win_sec: data_length-win_sec/2;
% t_ind = t_ind./(60*60); %%% change second scale to hour


T = config.T;
t_ind = 1 :T(end)/length(eig_info{1}): T(end);
t_ind = t_ind./(60); 



color_rgb_1 = [13, 55, 122]./255;
color_rgb_2 = [146, 15, 158]./255;

switch method_GC
    case 'Proposed'
        color_plot = color_rgb_1;
    case 'PNAS'
        color_plot = color_rgb_2;
end


max_eig_val = zeros(m(1),m(2));
sum_eig_val = zeros(m(1),m(2));

count = 0;
for f_ind = f_l:f_u
    count = count + 1;
    
    for win_num = 1 : m(2)
        eig_val_temp = diag(eig_val{count,win_num});
        max_eig_val(count,win_num) = eig_val_temp(1);
        sum_eig_val(count,win_num) = sum(eig_val_temp);
    end
    
end


%%  Plotting maximum of eigenvalues
count = 0;
for f_ind = f_l : f_u
    count = count + 1;
    figure
%     figure('units','normalized','outerposition',[0 0 1 1]),
    plot(t_ind , max_eig_val(count,:) , 'LineWidth',5 ,  'color' , 'r')
    xlabel('Time (mins)'),ylabel('Largest \lambda')
    xlim([0 t_ind(end)])

%     ylim([0 20000])
    set(gca,'FontSize',20)
    
%     str_tit_1 = sprintf('Maximum of Eigenvalues');
%     str_tit_2 = sprintf('%d channels at freq = %d *** window length = %d sec ** number of segment per each window = %d'...
%         ,ch_num, f_l ,win_sec , seg_num);
%     str_tit_3 = sprintf('%s' , method_GC);
%     title({str_tit_1 , str_tit_2 , str_tit_3});

end

%% Plotting summation of eigenvalues
count = 0;
for f_ind = f_l : f_u
     count = count + 1;  
    figure('units','normalized','outerposition',[0 0 1 1]),
    plot(t_ind , 20*log(sum_eig_val(count,:)) , 'LineWidth',5 ,  'color' , color_plot)
    xlabel('Time (hrs)'),ylabel('SUM of eigenvalues')
    xlim([0 t_ind(end)])

    str_tit_1 = sprintf('Summation of Eigenvalues');
    str_tit_2 = sprintf('%d channels at freq = %d *** window length = %d sec ** number of segment per each window = %d'...
        ,ch_num, f_l ,win_sec , seg_num);
    str_tit_3 = sprintf('%s' , method_GC);
    title({str_tit_1 , str_tit_2 , str_tit_3});

end


end