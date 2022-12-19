clc
clear
close all

Fs = 300;
Ts = 1/Fs;

% str_data = '8_12_2min_11973';  % not good 
% str_data = 'ali_reading_6_16_8529';
% str_data = 'ali_reading_6_16_eye_8526';
% str_data = 'reza_3';
% str_data = 'ali_0001';
% str_data = 'ali_11_29_22';
str_data = 'ali_11_30_22';

% str_data = '8_12_2min_4113';
% str_data = '8_12_2min_4699'; % repeating with above
% str_data = '8_12_2min_4707'; % repeating with above 

% str_data = 'A1'; % noise 
% str_data = 'A2'; % noise 

% str_data = 'ali_reading_6_16_[8529]'; % bad
% str_data = 'ali_reading_6_16_eye_[8526]'; 

str_load = sprintf('%s_raw.edf', str_data);
[hdr, data_raw] = edfread(str_load);

m = size(data_raw);

ch_num = m(1);
l_data = m(2);

t_sec = 0: Ts: l_data*Ts-Ts;
t_min = t_sec./60;

% rem_ch = [2, 3, 6, 11, 12, 19, 20, 22, 23, 24, 12, 13, 15, 16, 17, 20,   9, 17, 18, 21, 25, 26];
% rem_ch = [2, 3, 6, 19, 22, 24,    9, 17, 18, 21, 25, 26];
% rem_ch = [2, 3, 6, 22, 24,    9, 17, 18, 21, 25, 26];
rem_ch = [19, 9, 17, 18, 21, 25, 26];
% rem_ch = [0];
% rem_ch = [1, 2, 6, 7, 15, 16, 13, 10, 11, 12, 20, 23,   9, 17, 18, 21, 25, 26];
% rem_ch = [1, 2, 6, 7, 15, 16, 13, 10, 11, 12, 20, 23,   9, 17, 18, 21, 25, 26];


% for 8_12_2min_4113
% rem_ch = [11, 19,   9, 17, 18, 21, 25, 26];

% for ali_reading_6_16_eye_[8526] 
% rem_ch = [1, 10, 11, 12, 20,    9, 17, 18, 21, 25, 26];

cut_off_fr = 50/Fs/2; 
order = 16;
h = fir1(order, cut_off_fr);

data_det = zeros(ch_num-length(rem_ch), l_data);
cnt_ = 0;
for i=1: ch_num
    data_temp = data_raw(i, :);
    
    data_1 = [data_temp 0];
    data_2 = [0 data_temp];
    data_diff = data_1 - data_2;
    data_tmpt = data_diff(1: end-1);
    
    data_temp_rem_mean = data_temp-mean(data_temp);
    %     data_temp_rem_mean = data_temp;
    
%     data_temp_rem_mean_filt = conv(data_temp_rem_mean, h, 'same');
    
%     data_temp_rem_mean_filt = lowpass(data_temp_rem_mean, 5, Fs);
    
    data_i = data_temp_rem_mean;
    
    bw = 7;
    tpr_num = 9;
    params.tapers = [bw tpr_num];
    params.Fs = Fs;
    params.fpass = [5 35];
    
    movingwin = [4 2];
    [S,t,f] = mtspecgramc( data_i.', movingwin, params );
    
    
    

    
    if ~ismember(i, rem_ch)
        cnt_ = cnt_+1;
        data_det(cnt_, :) = data_i;
        disp(int2str(i));
        
        
        
        figure()
        subplot(121)
        imagesc(t/60, f, S.')
        colormap jet
        colorbar()
        set(gca,'YDir','normal')
        
        
        subplot(122)
        plot(data_i)
        
        str_ = sprintf('%d', i);
        suptitle(str_)
    end
end

strct_switch_data = [];
strct_switch_data.data_det = data_det;
strct_switch_data.Fs = Fs;
strct_switch_data.Ts = Ts;
strct_switch_data.t_sec = t_sec;

str_save = sprintf('%s.mat', str_data);
save(str_save, 'strct_switch_data');


u=1;



