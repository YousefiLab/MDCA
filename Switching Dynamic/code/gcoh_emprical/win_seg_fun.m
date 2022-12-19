function out = win_seg_fun(data_ch , config)
% out = WIN_SEG_FUN(data_ch , config) returns windowed data with inserted
% configuration
%
% INPUT
% data_ch       :  this matrix includes all channels of data 
%
% config is a structure with 9 fields in it
% 1 - ch_num    :  number of channels
% 2 - Fs        :  sampling rate of data
% 3 - sample_r  :  chosen sample rate for making windows and segments (this is not sample rate of data)
% 4 - win_sec   :  this value multiply to sample_r = length of window uses for chunking data
% 5 - win_length:  (sample) length of window that we want to calc Global Coherence on them
% 6 - seg_num   :  (number) number of segments that we want to devide a window into them  (number)
% 7 - seg_length:  (sample) number of samples of segments
% 8 - f_l       :  lower bound of desired freq interval
% 9 - f_u       :  upper bound of desired freq interval
% 10 - method_GC  specify the method that we use
% with method_GC, the method that we want use to calculate GC is selected
% there are two options,
% first = Proposed, means without subtracting from mean FFT
% second = PNAS, based PNAS Brown paper and subtracting FFT results from FFT's mean
% 11 - over_lap : (of 1) length of overlap window based win length 

%% Extract from config
ch_num = config.ch_num;
win_sec = config.win_sec;
win_length = config.win_length;
seg_num = config.seg_num;
sample_r = config.sample_r;
over_lap = config.over_lap;

signal_sample = length(data_ch); % length of signal


Fs = sample_r;
win_sample = win_sec*Fs;

%% calculate overlap indices
overlap_sample = floor(win_sample .* over_lap);

jump_sample = win_sample - overlap_sample;
 
% getting first and last indices of each window 
vec_step = 1:win_sample;
i = 0;
ex_con = 0;
while (ex_con == 0)
   i = i+1;
    temp_ind = (i-1)*jump_sample + vec_step;
    ind_end = temp_ind(end);
    if ind_end > signal_sample
        ex_con =1;
    end
    temp_ind_overlap(1,i) = temp_ind(1);
    temp_ind_overlap(2,i) = temp_ind(end);
    clear temp_indnum
end

ind_overlap = temp_ind_overlap(:,1:end-1);
% ind_overlap has 2 ROWS 
% first row indicates lower index of each overlapped window
% second row indicates upper index of each overlapped window
% COLUMNS shows value of first and last index number for different windows 
%% obtain data values of each window and segment
num_win = length(ind_overlap);


seg_sample = win_length./seg_num;
data_win_chunk = cell(1,ch_num);

for ch_num =1:ch_num

    %%% data_windowed = 1st index = 8s window number, 2nd index = times series of windowed data
    data_windowed = zeros(num_win,win_length);
    
    for i=1: num_win
        ind_l = ind_overlap(1,i);
        ind_u = ind_overlap(2,i);
        
        data_windowed(i,:) = data_ch(ind_l : ind_u , ch_num);
    end

    
    %%% data_win_chunk = 1st index = 8s window number, 2nd index = 1s chunk number, 3rd = time seires of windowed and chunk data
    data_win_seg_temp = zeros(num_win,seg_num,win_length./seg_num);
    
    for i=1:num_win
        data_win_temp = data_windowed(i,:);
        
        for j=1:seg_num
            data_win_seg_temp(i,j,:) = data_win_temp((j-1)*seg_sample + 1 : j*seg_sample);
        end
        
    end
    data_win_chunk{1,ch_num} = data_win_seg_temp;
end

out = data_win_chunk;
end