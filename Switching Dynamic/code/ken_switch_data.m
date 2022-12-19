clc
clear
close all

str_data = 'ken';
str_load = sprintf('%s.dat', str_data);
data = load(str_load);
data(:, 9) = [];
l_data = length(data);
ch_num = size(data, 2);

Fs = 300;
Ts = 1/Fs;

t_sec = 0: Ts: l_data*Ts-Ts;
t_min = t_sec./60;

ch_chosen = 1: 20;

% ch_chosen = [2, 5, 10, 18, 20];
data_det = zeros(length(ch_chosen), l_data);
cnt_ = 0;

for i=1: ch_num
    if ismember(i, ch_chosen)
        cnt_ = cnt_ + 1; 
        %    data_i = data(:, i).';
        %
        %    figure
        %    plot(t_min, data_i);
        %    xlabel('Time (min)');
        %    str_tit = sprintf('%d', i);
        %
        
        %
        %     data_temp = data(i, :);
        %
        %     data_1 = [data_temp 0];
        %     data_2 = [0 data_temp];
        %     data_diff = data_1 - data_2;
        %     data_tmpt = data_diff(1: end-1);
        %
        %     data_temp_rem_mean = data_temp-mean(data_temp);
        %     %     data_temp_rem_mean = data_temp;
        %
        %     %     data_temp_rem_mean_filt = conv(data_temp_rem_mean, h, 'same');
        %
        %     %     data_temp_rem_mean_filt = lowpass(data_temp_rem_mean, 5, Fs);
        %
        %     data_i = data_temp_rem_mean;
        
        
        data_i = data(:, i).';
        mean_i = mean(data_i);
        
        
        bw = 7;
        tpr_num = 9;
        params.tapers = [bw tpr_num];
        params.Fs = Fs;
        params.fpass = [5 35];
        
        movingwin = [4 2];
        [S,t,f] = mtspecgramc( data_i.', movingwin, params );
        
        
        data_det(cnt_, :) = data_i;
        disp(int2str(i));
        
        
        plt_= 0;
        if plt_ ==1
            figure()
            subplot(121)
            imagesc(t/60, f, S.')
            colormap jet
            colorbar()
            set(gca,'YDir','normal')
            
            
            subplot(122)
            plot(t_min, data_i)
            %     ind_ = find(ch_chosen==i);
            str_ = sprintf('%d', i);
            suptitle(str_)
        end
    end
    
end

strct_switch_data = [];
strct_switch_data.data_det = data_det;
strct_switch_data.Fs = Fs;
strct_switch_data.Ts = Ts;
strct_switch_data.t_sec = t_sec;
strct_switch_data.lbl_chosen = ch_chosen;

str_save = sprintf('%s.mat', str_data);
save(str_save, 'strct_switch_data');

u=1;