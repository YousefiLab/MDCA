function out = removed_fft_mean_fun(fft_win_chunk , config)

%%% config is a structure with 9 fields in it
%%% 1 - ch_num :    number of channels
%%% 2 - Fs :        sampling rate of data
%%% 3 - sample_r :  chosen sample rate for making windows and segments (this is not sample rate of data)
%%% 4 - win_sec :   this value multiply to sample_r = length of window uses for chunking data
%%% 5 - win_length  (sample) length of window that we want to calc Global Coherence on them
%%% 6 - seg_num     (number) number of segments that we want to devide a window into them  (number)
%%% 7 - seg_length  (sample) number of samples of segments
%%% 8 - f_l         lower bound of desired freq interval
%%% 9 - f_u         upper bound of desired freq interval
% 10 - method_GC  specify the method that we use
% with method_GC, the method that we want use to calculate GC is selected
% there are two options,
% first = Proposed, means without subtracting from mean FFT
% second = PNAS, based PNAS Brown paper and subtracting FFT results from FFT's mean

%%% Extract info
f_l = config.f_l;
f_u = config.f_u;
method_GC = config.method_GC;


m = size(fft_win_chunk);

fft_chunk_mean_removed = cell(1,m(2));
for ch_num = 1 : m(2)
    fft_chunk = cell2mat(fft_win_chunk(1,ch_num));
    
    mm=size(fft_chunk);
    fft_chunk_mean_removed_temp = zeros(mm(1) , mm(2) , f_u-f_l +1);
    count = 0;
    for f_index = f_l:f_u
        count = count + 1;
        for i=1:mm(1)
            
            chunk_temp = fft_chunk(i,:,f_index + 1);
            
            %%% different methods
            switch method_GC
                
                case 'Proposed'
                    mean_fft = 0;
                case 'PNAS'
                    mean_fft = mean(chunk_temp);
            end
            
            fft_chunk_mean_removed_temp(i,:,count) = chunk_temp -mean_fft;
            
            
        end
    end
    fft_chunk_mean_removed{1,ch_num} = fft_chunk_mean_removed_temp;
end

out = fft_chunk_mean_removed;
end