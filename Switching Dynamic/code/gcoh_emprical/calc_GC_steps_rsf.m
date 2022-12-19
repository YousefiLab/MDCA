function [GC_overal , sorted_eig_info, fft_seg_mean_removed] = calc_GC_steps_rsf(data_chs , config)
% CALC_GC_STEPS calculate golbal coehrence based on input configuration and
% after calculating GC, the value of it, decomposed to EIGEN componnets and
% then sort EIGENVALUEs descendingly and after that sorts each EIGENVECTOR
% corrosponded to its eigenvalue as a column of eigenvector martrix
%
%Input 
% 1- data_chs a matrix indicates input time series
% columns = different channels
% rows = samples of each channel
%
%2-config
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
% 12 - method_angle :choosing method of calculating angle for complex vectors, there is 2 options
%
% OUTPUT
% 1- GC : a matrix that returns GC value for desired frequencies and desired
% window numbers
% 2- sorted_eig_info : a cell with size 2 
% first cell = matrix of EIGENVALUES after sorting descendingly
% second cell = matrix of EIGENVECTORS after sorting descendingly

%% STEP 1 (WINDOWING DATA)

%%% this fun devides data into desired windows and then each window devides
%%% into desired segments
data_win_seg = win_seg_rsf(data_chs , config);
 
%% STEP 2 (TAPERED FFT)

%%% this part of code produce a hanning window and multiply it into segments
%%% of data and after that calculate fft
fft_win_seg = slepian_fft_rsf(data_win_seg , config);

%% STEP 3 (CROSS SPECTRAL MATRIX)

%%% this fucntion calculate fft at desired freq interval and then subtract
%%% calculated fft from mean of diffrenet segments

fft_seg_mean_removed = removed_fft_mean_rsf(fft_win_seg , config);


%%% this fun calculate cross spectra at desired freqs for all windows
%%% the result is a cell,
%%% number of rows show number of desired freq
%%% number of columns show number of windows


cross_spect = cross_spect_rsf(fft_seg_mean_removed);

%% STEP 4 (eigenvalueS AND VECTORS and CALC GLOBAL COHERENCE)

%%% first output = a cell that includes 2 cell;
%%% first of those =  stored eigenvalue's matrix of different windows (col) and diffrent frequncies(row)
%%% second = stored eigen vector's matrix of different windows (col) and diffrent frequncies(row)

%%% second output = a matrice that shows Global Coherence value for different windows (col) and different freqs (row)
% % [Eig_info , GC_slep, GC_overal] = eig_rsf(cross_spect);
[eig_info , GC_overal] = eig_fun(cross_spect);

%%% this fun sorting eigvals and corresponding eigvecs in a descend way
% sorted_eig_info = sort_eig_info_rsf(Eig_info);
sorted_eig_info = sort_eig_info_fun(eig_info);


%% SAVE INFO 

str_save = sprintf('anesthesia_ch%d_win%d_overLap%.0f_fL%d_fU%d_numSlep%d.mat',...
    config.ch_num , config.win_sec, 100*config.over_lap, config.f_l, config.f_u, config.num_slepian);
save(str_save , 'GC_overal', 'sorted_eig_info', 'fft_seg_mean_removed');
end