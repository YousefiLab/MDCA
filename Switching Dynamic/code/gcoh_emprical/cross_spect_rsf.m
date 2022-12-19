function out = cross_spect_rsf(FFT_REM, win_gc, over_lap_gc)


% length of jumping with considering overlap
% j_s = win_gc - over_lap_gc;  % Jumping Step


% number of slepian
pp = size(FFT_REM);

num_ch = size(FFT_REM{1});              %% number of channels

%%% m(1) = win's number ** m(2) = segment's number
%%% m(3) = desired freq number (if there is not m(3) we only have one desired freq)
m = size(cell2mat(FFT_REM{1}(1)));



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
    for win_num=1: m(1)
        
        %%% calc CROSS SPECTRAL between chaanals
        cross_spect_temp = zeros(num_ch(2),num_ch(2));
        for i=1: num_ch(2)
            vec_fft_temp = [];
            for slep_num=1: pp(2)
                FFT_REM_temp = FFT_REM{slep_num};
                vec_fft_temp = [vec_fft_temp FFT_REM_temp{i}(win_num , : , f_ind)];
            end
                vec_fft(i,:) = vec_fft_temp;
        end
        
        
        
        %% calculate cross spectral matrix
        
        for kk=1:num_ch(2)
            for ll=1:num_ch(2)
                aa = vec_fft(kk,:);
                bb = vec_fft(ll,:)';
                cross_spct_temp(kk,ll) = (vec_fft(kk,:)*vec_fft(ll,:)');
            end
        end
        
        
        
        cross_spect{f_ind, win_num} = cross_spct_temp;
        
    end
end

out = cross_spect;


end