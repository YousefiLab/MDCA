function out = angle_eig_vec_fun(eig_info , config, fnt_size)
%this function calculate angle between largest correspondin eigenvalues for
%2 steps k and k+1
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
% 11 - over_lap : (of 1) length of overlap window based win length
% 12 - method_angle :choosing method of calculating angle for complex vectors, there is 2 options
% first = 1  : uses real part of complex vectors
% second = 2 : uses Hermitain angle
% there is more information for these method in following link:
% https://www.researchgate.net/profile/Wiwat_Wanicharpichat/post/How_can_I_calculate_the_angle_between_two_complex_vectors/attachment/59d645cfc49f478072eae128/AS%3A273828859580416%401442297296461/download/9904077.pdf
%
%
%OUTPUT
% The output is a vector with size of number of window -1, that returns the
% Hermitain angle between step k and k+1 largest crresponding eigenvectors

%%%% Extract config info
ch_num = config.ch_num;
sample_r = config.sample_r;
win_length = config.win_length;
seg_num = config.seg_num;
f_l = config.f_l;
method_GC = config.method_GC;
method_angle = config.method_angle;

eig_vec = eig_info{2};

m = size(eig_vec);
des_fre_num = m(1);
win_num = m(2);
%%% des_fre_num number of desired freqs
%%% win_num number of windows

win_sec = win_length./sample_r;
% data_length = win_num*win_sec;
% t_ind = win_sec/2 :win_sec: data_length-win_sec/2;
% t_ind = t_ind./(60*60); %%% change second scale to hour


T = config.T;
t_ind = 1 :T(end)/length(eig_info{1}): T(end);
t_ind = t_ind./(60); 


% because the length of angle between eigenvalues is one less than number
% of windows, the first inex is elminated
t_ind = t_ind(2:end);

color_rgb_1 = [174, 16, 232]./255;
color_rgb_2 = [244, 209, 66]./255;



switch method_GC
    case 'Proposed'
        color_plot = color_rgb_1;
    case 'PNAS'
        color_plot = color_rgb_2;
end

%% relying REAL vectors
if method_angle == 1
    
    THETA_REAL = zeros(des_fre_num , win_num-1);
    for fr_num = 1 : des_fre_num
        
        for k = 1: win_num-1
            
            %%% columns are corresponding eigevectors
            temp_eig_vec_i = eig_vec{fr_num , k}(:,1);     %% largest corresponding eigenvector for step i
            real_temp_eig_vec_i = real(temp_eig_vec_i);
            imag_temp_eig_vec_i = imag(temp_eig_vec_i);
            
            temp_eig_vec_j = eig_vec{fr_num , k+1}(:,1);  %% largest corresponding eigenvector for step i+1
            real_temp_eig_vec_j = real(temp_eig_vec_j);
            imag_temp_eig_vec_j = imag(temp_eig_vec_j);
            
            length_vec = length(real_temp_eig_vec_j);
            relying_real_eig_vec_i = zeros(2*length_vec , 1);
            relying_real_eig_vec_j = zeros(2*length_vec , 1);
            
            for i=1 : length_vec
                relying_real_eig_vec_i(2*i-1 , 1) = real_temp_eig_vec_i(i);
                relying_real_eig_vec_i(2*i , 1) = imag_temp_eig_vec_i(i);
                
                relying_real_eig_vec_j(2*i-1 , 1) = real_temp_eig_vec_j(i);
                relying_real_eig_vec_j(2*i , 1) = imag_temp_eig_vec_j(i);
            end
            
            
            CosTheta = dot(relying_real_eig_vec_i , relying_real_eig_vec_j)/(norm(relying_real_eig_vec_i)*norm(relying_real_eig_vec_j));
            
            THETA_REAL(fr_num , k) = acos(CosTheta);
            
        end
        
    end
    out = THETA_REAL;
end
%% Hermitian Method
if method_angle == 2
    
    THETA_HERMITIAN = zeros(des_fre_num , win_num-1);
    for fr_num = 1 : des_fre_num
        
        for k = 1: win_num-1
            
            %%% columns are corresponding eigevectors
            temp_eig_vec_i = eig_vec{fr_num , k}(:,1);     %% largest corresponding eigenvector for step i
            temp_eig_vec_j = eig_vec{fr_num , k+1}(:,1);  %% largest corresponding eigenvector for step i+1
            
            CosTheta = dot(temp_eig_vec_i , temp_eig_vec_j)/(norm(temp_eig_vec_i)*norm(temp_eig_vec_j));
            
            psudo_angle = angle(CosTheta); %%% -pi <=psudo_angle <=pi
            ro_hermitian = abs(CosTheta);  %%% 0 <= ro_hermitian <= 1
            
            THETA_HERMITIAN(fr_num , k) = acos(ro_hermitian);
            
        end
        
    end
    
    out = THETA_HERMITIAN;
    
end


%% PLOTTING PART
switch method_angle
    case 1
        str_method = 'Relying Real Angle (rad)';
        
    case 2
        str_method = 'Hermitan Angle (rad)';
end

for fr_num =1 : des_fre_num
    temp_theta_eig_vec = out(fr_num , :);
    
%     figure('units','normalized','outerposition',[0 0 1 1]),
    
    figure('units','normalized','outerposition',[0 0 1 1]);

    
    
%     theta_plot = abs(temp_theta_eig_vec);
    theta_plot = (temp_theta_eig_vec);
    plot(t_ind , theta_plot , 'LineWidth', 2 , 'color' , 'b')
    xlim([0 t_ind(end)]) , 
    ylim([0 pi./2]);
    xlabel('Time (mins)');
    ylabel({'Hermitian'; 'Angle (rad)'});
    
    set(gca,'FontSize', fnt_size)
    
    set(gca, 'YTIck', [0:.5:1.9])
    set(gca, 'YAxisLocation', 'right')
    set(gca, 'XTIck', [20:40:140])
    
    set(gcf, 'PaperPosition', [0 0 6 4]);
    print('herm_Larg_12','-dpng','-r600')
    
    
%     str_tit_1 = sprintf('Angle between Largest Corresponding Eigenvectors');
%     str_tit_2 = sprintf('%d channels at freq = %d *** window length = %d sec ** number of segment per each window = %d'...
%         ,ch_num, f_l ,win_sec , seg_num);
%     str_tit_3 = sprintf('%s' , method_GC);
%     title({str_tit_1 , str_tit_2 , str_tit_3});
    
end


end