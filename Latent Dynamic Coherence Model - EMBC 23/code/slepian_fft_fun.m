function out = slepian_fft_fun(data_win_seg , config)

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

%%% Extract Info
seg_length = config.seg_length;
sample_r = config.sample_r;

m = size(data_win_seg);

%%% slapeian filter
seq_length = floor(seg_length);
time_halfbandwidth = 2.5;
num_seq = 1;
[dps_seq] = dpss(seq_length,time_halfbandwidth,num_seq);
dps_seq = dps_seq.'./sqrt(sum(dps_seq.^2));

% figure
% plot(dps_seq, 'LineWidth' , 2 , 'color' , [62, 140, 168]./255)
% xlim([1 128])
% xlabel('sample')
% title('Slepian Sequences, N = 128, NW = 2.5' );grid on


slepian_taper = dps_seq;


fft_win_seg = cell(1,m(2));
for ch_num=1 : m(2)
    data_win_seg_temp = cell2mat(data_win_seg(1,ch_num));
    mm = size(data_win_seg_temp);
    fft_length = sample_r;
    fft_win_seg_temp = zeros(mm(1) , mm(2) , mm(3));
    for i=1: mm(1)
        for j=1:mm(2)
            data_temp = data_win_seg_temp(i,j,:);
            data_temp = reshape(data_temp , 1,floor(seg_length)); %%% convert from 3d mat into 2d mat
            slep_data = slepian_taper.*data_temp;
            
            
            
            Y = fft(slep_data , fft_length);
            fft_win_seg_temp(i,j,:) = Y(1:mm(3));
            
%             disp(j)
            figure
            plot(abs(Y(1:mm(3))), 'LineWidth' , 3 , 'color' , [90, 108, 230]./255)
            j
%             xlim([1 128])
%             xlabel('freq (Hz)')
%             title('FFT of tapered segment 16' )
%             close all
            
            
            %%%%%%%% Test
% %             Y = Y(1:mm(3));
% %             
% %             figure,
% %             stem(abs((Y)))
% %             close
            %%%%%%%%%% test
        end
    end
    fft_win_seg{1,ch_num} = fft_win_seg_temp;
end

out = fft_win_seg;
end