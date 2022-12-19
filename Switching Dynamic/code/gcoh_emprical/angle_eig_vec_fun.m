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
method_angle = config.method_angle;

eig_vec = eig_info{2};

m = size(eig_vec);
des_fre_num = m(1);
win_num = m(2);
line_width = 2;



T = config.t_min;
step_t = T(end)/length(eig_info{1});
t_ind = 0 :step_t :T(end)-step_t;


% because the length of angle between eigenvalues is one less than number
% of windows, the first inex is elminated
t_ind = t_ind(2:end);

color_rgb_1 = [174, 16, 232]./255;
color_rgb_2 = [244, 209, 66]./255;

%% calc angle
switch method_angle 
    case 1 % relying REAL vectors
        out_largest= calc_angle_meth_1(eig_vec, des_fre_num, win_num, 1);
        out_sec_largest= calc_angle_meth_1(eig_vec, des_fre_num, win_num, 2);
    case 2 % Hermitian Method
        out_largest = calc_angle_meth_2(eig_vec, des_fre_num, win_num, 1);
        out_sec_largest = calc_angle_meth_2(eig_vec, des_fre_num, win_num, 2);
end

%% PLOTTING PART
for fr_num =1 : des_fre_num % for all freq
    temp_theta_eigVec_largest = out_largest(fr_num , :);
    temp_theta_eigVec_sec = out_sec_largest(fr_num , :);

    left_pos = .25;
    bot_pos = .25;
    width_val = .7;
    height_val = .33;
    
    bot_diff = .37;
    
    h=figure;
    %% second largest
    ax = subplot(212);
    ax.Position = [left_pos (bot_pos+ 0*bot_diff) width_val height_val];
    
    theta_plot = (temp_theta_eigVec_sec);
    plot(t_ind , theta_plot , 'LineWidth', line_width , 'color' , 'r')
    ylim([0 pi./2]);
    yticks([0 pi./4])
    yticklabels({'0','\pi/4'})

    
    xlabel('Time (mins)');
%     ylabel({'Hermit Ang (rad) '});
    
    set(gca,'FontSize', fnt_size)
    
    % xtick
    ax =gca;
    x_num_step = 3;
    x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
    ax.XTick = x_tick_values;

%     set(gcf, 'PaperPosition', [0 0 6 4]);
%     print('herm_Larg_12','-dpng','-r600')
    

   %% Second Largest Plot 

    ax = subplot(211);
    ax.Position = [left_pos (bot_pos+ 1*bot_diff) width_val height_val];
    
    theta_plot = (temp_theta_eigVec_largest);
    plot(t_ind , theta_plot , 'LineWidth', line_width , 'color' , 'b')
    xlim([0 t_ind(end)]), 
    ylim([0 pi./2]);
    yticks([0 pi./4])
    yticklabels({'0','\pi/4'})
    
    ax = gca;
    ax.TickLabelInterpreter
    xticklabels({'0','\pi\2'})
    %     xlabel('Time (mins)');
    %     ylabel({'Hermit Ang (rad) '});
    hh = text(-30, -1.5,'Hermitian Angle (rad)', 'FontSize', fnt_size);
    set(hh,'Rotation',90);
    
    set(gca,'FontSize', fnt_size)
    
    % xtick
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
%     ax =gca;
%     x_num_step = 3;
%     x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
%     ax.XTick = x_tick_values;

%     set(gcf, 'PaperPosition', [0 0 6 4]);
%     print('herm_Larg_12','-dpng','-r600')

    
end


end


%% nested functions 
% method 1 - relying REAL vectors
function out= calc_angle_meth_1(eig_vec, des_fre_num, win_num, eig_num)
    
THETA_REAL = zeros(des_fre_num , win_num-1);
for fr_num = 1 : des_fre_num
    
    for k = 1: win_num-1
        
        %%% columns are corresponding eigevectors
        temp_eig_vec_i = eig_vec{fr_num , k}(:,eig_num);     %% largest corresponding eigenvector for step i
        real_temp_eig_vec_i = real(temp_eig_vec_i);
        imag_temp_eig_vec_i = imag(temp_eig_vec_i);
        
        temp_eig_vec_j = eig_vec{fr_num , k+1}(:,eig_num);  %% largest corresponding eigenvector for step i+1
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

% method 2- hermitian angle
function out = calc_angle_meth_2(eig_vec, des_fre_num, win_num, eig_num)

THETA_HERMITIAN = zeros(des_fre_num , win_num-1);
for fr_num = 1 : des_fre_num
    
    for k = 1: win_num-1
        
        %%% columns are corresponding eigevectors
        temp_eig_vec_i = eig_vec{fr_num , k}(:,eig_num);     %% largest corresponding eigenvector for step i
        temp_eig_vec_j = eig_vec{fr_num , k+1}(:,eig_num);  %% largest corresponding eigenvector for step i+1
        
        CosTheta = dot(temp_eig_vec_i , temp_eig_vec_j)/(norm(temp_eig_vec_i)*norm(temp_eig_vec_j));
        
        psudo_angle = angle(CosTheta); %%% -pi <=psudo_angle <=pi
        ro_hermitian = abs(CosTheta);  %%% 0 <= ro_hermitian <= 1
        
        THETA_HERMITIAN(fr_num , k) = acos(ro_hermitian);
        
    end
    
end

out = THETA_HERMITIAN;


end



% finding the proper val for ticks
function tick_val = fun_finding_tick(ax, num_step, axis)

switch axis
    case 'x'
        lim_val = ax.XLim;
    case 'y'
        lim_val = ax.YLim;
end

step_tick = ((lim_val(end) - lim_val(1))./num_step);
cnt = find_order(step_tick); % nested fun
if cnt> 1
    step_tick_n = round(step_tick, -(cnt-1));
    tick_val = round(lim_val(1) + (step_tick_n./2), -(cnt-1)) : step_tick_n: lim_val(end);
else
    step_tick_n = round(step_tick, cnt);
    tick_val = round(lim_val(1) + (step_tick_n./2), cnt) : step_tick_n: lim_val(end);
end


end

% finding the order of input value
function cnt = find_order(x)

if x> 1
    stp_con = 0;
    cnt = -1;
    while stp_con == 0
        cnt = cnt +1;
        a = floor(x./10^cnt);
        
        if a == 0
            stp_con=1;
        end
        
    end
else
    stp_con = 0;
    cnt = 0;
    while stp_con == 0
        cnt = cnt +1;
        a = floor(x./10^(-cnt));
        
        if a > 0
            stp_con=1;
        end
        
    end
end

end
