function out = cross_spect_fun(FFT_REM)

fft_seg_mean_removed_ch = cell2mat(FFT_REM(1));

num_ch = size(FFT_REM);              %% number of channels

%%% m(1) = win's number ** m(2) = segment's number 
%%% m(3) = desired freq number (if there is not m(3) we only have one desired freq)
m = size(fft_seg_mean_removed_ch);   



%%% calculate number of desired freqs
if length(m)>2 
    ind_end = m(3);
else
    ind_end = 1;
end

%%% cross_spect cell indicates cross spectral of different freqs and different windows
%%% the value of rows shows different desired freqs and columns shows different time windows
cross_spect = cell(ind_end , m(1));


for f_ind = 1:ind_end

    %%% interval of desired freqs
    for win_num = 1: m(1)

        %%% calc CROSS SPECTRAL between chaanals
        cross_spect_temp = zeros(num_ch(2),num_ch(2));
        for i=1:num_ch(2)
            for j=1:num_ch(2)
                
                X_i = FFT_REM{i}(win_num,:,f_ind);
                X_j = FFT_REM{j}(win_num,:,f_ind);
                K = length(X_i);

               cross_spect_temp(i,j) = (1/K)*sum(X_i.*conj(X_j));
            end
        end
        
        cross_spect{f_ind,win_num} = cross_spect_temp;
        
    end
end

out = cross_spect;
end