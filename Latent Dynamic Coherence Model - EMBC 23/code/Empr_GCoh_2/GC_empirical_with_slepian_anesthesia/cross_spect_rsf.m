function out = cross_spect_rsf(FFT_REM, win_gc, over_lap_gc)


% length of jumping with considering overlap
% j_s = win_gc - over_lap_gc;  % Jumping Step


% number of slepian
pp = size(FFT_REM);
Num_slep = pp(2);

num_ch = size(FFT_REM{1}, 2);              %% number of channels

%%% m(1) = win's number ** m(2) = segment's number
%%% m(3) = desired freq number (if there is not m(3) we only have one desired freq)
m = size(cell2mat(FFT_REM{1}(1)));
num_chunk = m(1);


%%% calculate number of desired freqs
if length(m)>2
    ind_end = m(3);
else
    ind_end = 1;
end

%%% cross_spect cell indicates cross spectral of different freqs and different windows
%%% the value of rows shows different desired freqs and columns shows different time windows
cross_spect = cell(ind_end , num_chunk);


for f_ind = 1:ind_end
    
    
    %%% interval of desired freqs
    for win_num=1: num_chunk
        
        %%% calc CROSS SPECTRAL between chaanals
        cross_spect_temp = zeros(num_ch,num_ch);
        for ch_=1: num_ch
            vec_fft_temp = [];
            for tpr=1: Num_slep
                FFT_REM_temp = FFT_REM{tpr};
                a = FFT_REM_temp{ch_}(win_num , : , f_ind);
                vec_fft_temp = [vec_fft_temp a];
            end
                vec_fft(ch_,:) = vec_fft_temp;
        end
        
        
        
        %% calculate cross spectral matrix
        
        for kk=1:num_ch
            for ll=1:num_ch
                aa = vec_fft(kk,:);
                bb = vec_fft(ll,:)';
                cc = (vec_fft(kk,:)*vec_fft(ll,:)');
                cross_spct_temp(kk,ll) = cc;
            end
        end
        
        
        
        cross_spect{f_ind, win_num} = cross_spct_temp;
        
    end
end

out = cross_spect;


end