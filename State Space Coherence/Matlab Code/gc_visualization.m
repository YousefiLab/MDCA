function [] = gc_visualization(param_v)
close all
clear
clc

load('L_update_iter15.mat')

REAL_PARAM = param_v.REAL_PARAM;
PARAM = param_v.PARAM;
ParamaUpdate = param_v.ParamaUpdate;
pdf_x = param_v.pdf_x;
x = param_v.x;
x_k = param_v.x_k;

% number of iteration
iter_num = length(PARAM);

% get RMSE and the value of mean and max x of distribution
[X_mean, X_max, rmse_mean, rmse_max] = calculate_rmse(pdf_x, x, x_k);

sv_update = ParamaUpdate.sv_update; 
mu_update = ParamaUpdate.mu_update;
L_update = ParamaUpdate.L_update;
ab_update = ParamaUpdate.ab_update;


% s_ab_real = REAL_PARAM.s_ab;
sv_real = REAL_PARAM.sv;
mu_real = REAL_PARAM.mu;
L_real = REAL_PARAM.L;
am_real = REAL_PARAM.ab.a;
bm_real = REAL_PARAM.ab.b;

iter_num = length(PARAM);

% obtain data of updated parameters
for iter=1: iter_num
    
    if sv_update ==1 
        sv_estimate(1,iter) = PARAM{iter}.sv;
    end
    if mu_update==1
        mu_estimate(:,iter) = PARAM{iter}.mu;
    end
    
    if L_update==1
        L_estimate(:, :, iter) = PARAM{iter}.L;
    end
    
    if ab_update
        am_estimate(:,iter) = PARAM{iter}.ab(:,1);
        bm_estimate(:,iter) = PARAM{iter}.ab(:,2);
    end
    
end


% cal plot functions
for iter=1: iter_num-1
    if sv_update ==1 
        plot_sv(sv_estimate, sv_real, iter_num)
    end
    
    count_update_all = sv_update+ mu_update+ L_update + ab_update;
    %% only updated mu 
    if mu_update ==1 && count_update_all==1
        
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        
        % rmse
        color_1 = [204, 0, 82]./255;
        color_2 = [31, 122, 122]./255;
        pv_rmse_1 = [0.02, 0.60, 0.2, 0.35];
        pv_rmse_2 = [0.24, 0.60, 0.27, 0.35];
        plot_rmse(rmse_mean, rmse_max, X_mean{iter}, X_max{iter}, x_k, h, iter, pv_rmse_1, pv_rmse_2, 10, color_1, color_2)
        
        
        % pdf
        pv_pdf = [0.02, 0.1, 0.49, 0.4];
        plot_pdf_x(pdf_x, x, x_k, iter, pv_pdf,h)
        
        % mu
        color_real = [255, 77, 166]./255;
        color_imag = [204, 51, 0]./255;
        pv_mu1 = [.56 ,.6 ,.2, .35];
        pv_mu2 = [.79 ,.6 ,.2, .35];
        pv_mu3 = [.56 ,.1 ,.2, .35];
        pv_mu4 = [.79 ,.1 ,.2, .35];
        plot_mu(mu_estimate, mu_real, iter_num, pv_mu1, pv_mu2, pv_mu3, pv_mu4,h, iter, 10, color_real, color_imag)

        % save part
        str_save_mu_dist = sprintf('mu_rmse_distX_iter_%d.jpg', iter);
        saveas(h, str_save_mu_dist);
        close
    end
    
    if L_update==1 && count_update_all==1
                h = figure('units','normalized','outerposition',[0 0 1 1]);
        
        % rmse
        color_1 = [204, 0, 82]./255;
        color_2 = [31, 122, 122]./255;
        pv_rmse_1 = [0.02, 0.60, 0.2, 0.35];
        pv_rmse_2 = [0.24, 0.60, 0.27, 0.35];
        plot_rmse(rmse_mean, rmse_max, X_mean{iter}, X_max{iter}, x_k, h, iter, pv_rmse_1, pv_rmse_2, 10, color_1, color_2)
        
        
        % pdf
        pv_pdf = [0.02, 0.1, 0.49, 0.4];
        plot_pdf_x(pdf_x, x, x_k, iter, pv_pdf,h)
        
        % L 
        pv_real_L = [.6,.6,.3,.35];
        pv_imag_L = [.6,.1,.3,.35];
        plot_L(L_estimate, L_real, pv_real_L, pv_imag_L, h, iter )
        
        %         save part
        str_save_mu_dist = sprintf('L_rmse_distX_iter_%d.jpg', iter);
        saveas(h, str_save_mu_dist);
        close
    end
    
    %% only updated am, bm
    if ab_update==1 && count_update_all==1
        
        h = figure('units','normalized','outerposition',[0 0 1 1]);

        % rmse
        pv_rmse_1 = [0.02, 0.60, 0.2, 0.35];
        pv_rmse_2 = [0.24, 0.60, 0.27, 0.35];
        color_rmse_1 = [204, 0, 82]./255;
        color_rmse_2 = [31, 122, 122]./255;
        plot_rmse(rmse_mean, rmse_max, X_mean{iter}, X_max{iter}, x_k, h, iter, pv_rmse_1, pv_rmse_2, 10, color_rmse_1, color_rmse_2)
        
        % pdf
        pv_pdf = [0.02, 0.1, 0.49, 0.4];
        plot_pdf_x(pdf_x, x, x_k, iter, pv_pdf,h)
        
        % rmse
        pv_rmse_1 = [0.02, 0.60, 0.2, 0.35];
        pv_rmse_2 = [0.24, 0.60, 0.27, 0.35];
        color_rmse_1 = [204, 0, 82]./255;
        color_rmse_2 = [31, 122, 122]./255;
        plot_rmse(rmse_mean, rmse_max, X_mean{iter}, X_max{iter}, x_k, h, iter, pv_rmse_1, pv_rmse_2, 10, color_rmse_1, color_rmse_2)
        
        % am, bm
        pv_am = [.56 ,.1 ,.4, .35];
        pv_bm = [.56 ,.6 ,.4, .35];
        color_am = [0, 77, 0]./255;
        color_bm = [153, 51, 51]./255;
        plot_am_bm(am_estimate, bm_estimate, am_real, bm_real, iter_num, iter, h, pv_am, pv_bm, 10, color_am, color_bm)
        
        % save part
        str_save_ab_dist = sprintf('ab_rmse_distX_iter_%d.jpg', iter);
        saveas(h, str_save_ab_dist);
        close
    end
    
    %% simultanously updated am, bm and my
    if ab_update*mu_update ==1  && count_update_all==2
               
        h = figure('units','normalized','outerposition',[0 0 1 1]);
       
        
        % rmse
        pv_rmse_1 = [0.02, 0.60, 0.2, 0.35];
        pv_rmse_2 = [0.24, 0.60, 0.27, 0.35];
        color_rmse_1 = [204, 0, 82]./255;
        color_rmse_2 = [31, 122, 122]./255;
        plot_rmse(rmse_mean, rmse_max, X_mean{iter}, X_max{iter}, x_k, h, iter, pv_rmse_1, pv_rmse_2, 10, color_rmse_1, color_rmse_2)
        
        
        % am, bm
        pv_am = [0.04 , .1, .22, .35];
        pv_bm = [.3 , .1, .22, .35];
        color_am = [0, 77, 0]./255;
        color_bm = [153, 51, 51]./255;
        plot_am_bm(am_estimate, bm_estimate, am_real, bm_real, iter_num, iter, h, pv_am, pv_bm, 10, color_am, color_bm)
        
        
        % mu
        color_real = [255, 77, 166]./255;
        color_imag = [204, 51, 0]./255;
        pv_mu1 = [.56 ,.6 ,.2, .35];
        pv_mu2 = [.79 ,.6 ,.2, .35];
        pv_mu3 = [.56 ,.1 ,.2, .35];
        pv_mu4 = [.79 ,.1 ,.2, .35];
        plot_mu(mu_estimate, mu_real, iter_num, pv_mu1, pv_mu2, pv_mu3, pv_mu4,h, iter, 10, color_real, color_imag)
        
        % save part
        str_save_ab_dist = sprintf('ab_mu_rmse_distX_iter_%d.jpg', iter);
        saveas(h, str_save_ab_dist);
        close
%         
    end
    
    
    
end




end

%% RMSE CALCULATES
function [X_mean, X_max, rmse_mean, rmse_max] = calculate_rmse(pdf_x, x, x_k)

ITER = length(pdf_x);
K = length(x_k);
for i=1:ITER
    temp_pdf = pdf_x{i};

    
    x_mean = zeros(K,1);
    x_max = zeros(K,1);
    for k=1:K
        %     [pdf_sorted, I] = sort(temp_pdf(:,k));
        pdf_k = temp_pdf(:,k);
        mean_pdf = mean(pdf_k);
        max_pdf = max(pdf_k);
        
        ind_u_mean = find(pdf_k>=mean_pdf); % index Uper mean
        x_mean(k,1) = x(ind_u_mean(1));
        
        ind_max = find(pdf_k>=max_pdf); % index Uper mean
        x_max(k,1) = x(ind_max(1));
        
    end
    
    X_mean{i} = x_mean;
    X_max{i} = x_max;
    
    rmse_mean(i) = sqrt(sum((x_k-x_mean).^2)./K);
    rmse_max(i) = sqrt(sum((x_k-x_max).^2)./K);
    

end


end

%% plot
% RMSE
function [] = plot_rmse(rmse_mean, rmse_max, x_mean, x_max, x_k, h, iter, pv_1, pv_2, color_size, color_1, color_2)

iter_num = length(rmse_mean);
K = length(x_k);
%% 1
axes1 = axes('Parent',h,'Position',pv_1);
hold(axes1,'on');
ylim_min = min(min(rmse_mean),min(rmse_max))-.5;
ylim_max = max(max(rmse_mean), max(rmse_max))+.5;
ylim(axes1,[ylim_min ylim_max]);
xlim(axes1,[1 iter_num]);

plot(1:iter_num, rmse_mean, 'Color', color_1,'LineWidth', 1.5), hold on 
plot(1:iter_num, rmse_max, 'Color', color_2,'LineWidth', 1.5), hold on 
xlabel('iter')
legend('RMSE_{mean}', 'RMSE_{max}')

hold on
p = plot(axes1 , iter, rmse_mean(iter) , 'o');
p.MarkerSize = color_size;
p.MarkerFaceColor = color_1;
p.MarkerEdgeColor = color_1;

hold on
p = plot(axes1 , iter, rmse_max(iter) , 'o');
p.MarkerSize = color_size;
p.MarkerFaceColor = color_2;
p.MarkerEdgeColor = color_2;

%% 2
axes2 = axes('Parent',h,'Position',pv_2);
hold(axes2,'on');
% ylim_min = min([min(x_k) min(x_mean) min(x_max)])-.5;
% ylim_max = max([max(x_k) max(x_mean) max(x_max)])+.5;
ylim_min = min(x_k) -10;
ylim_max = max(x_k) +2;

ylim(axes2,[ylim_min ylim_max]);
xlim(axes2,[1 K]);

plot(1:K, x_k, 'k','LineWidth', 2.5), hold on
plot(1:K, x_mean, 'Color', color_1,'LineWidth', 1.5), hold on
plot(1:K, x_max, 'Color', color_2,'LineWidth', 1.5)
xlabel('step')
legend('x_k', 'x_{min}', 'x_{max}')



    
end

% pdf 
function [] = plot_pdf_x(pdf_x, x, x_k, iter, positionVector, h)

K = length(x_k);


    pdf_x_temp = pdf_x{iter};
    
%     subplot('Position', positionVector)
    
    axes1 = axes('Parent',h, 'Position',positionVector);
    hold(axes1,'on');
    
    imagesc(1:K, x, pdf_x_temp), hold on
    plot(1:K, x_k,'r', 'LineWidth',2);
    xlabel('step')
    ylim(axes1,[min(x) max(x)]);
    xlim(axes1,[1 K]);

    title(['estimated distribution of x and x_k - iter = ', num2str(iter)])
    



end

% plot am & bm
function [] = plot_am_bm(am_updated, bm_updated, am_real, bm_real, iter_num, iter, h, pv_am, pv_bm, color_size, color_am, color_bm)
s_am = size(am_updated);
% s_am(1)   : channel num
% s_am(2)   : iter

for i=2: s_am(1)
    
    %% am
    am_update_ch = am_updated(i,:);
    
    axes1 = axes('Parent',h,'Position',pv_am);
    hold(axes1,'on');
%     ylim(axes1,[min(data)-.5 max(data)+.5]);
    xlim(axes1,[1 iter_num]);
    
    plot(1:iter_num, am_update_ch, 'Color', color_am, 'LineWidth', 1.5),hold on
    plot(1:iter_num, am_real(i).*ones(1,iter_num), 'b', 'LineWidth', 1.5)
    title('a_m')
%     legend('estimated value', 'real value')
    
    hold on
    p = plot(axes1 , iter, am_update_ch(iter) , 'o');
    p.MarkerSize = color_size;
    p.MarkerFaceColor = color_am;
    p.MarkerEdgeColor = color_am;
    
    %% bm
    bm_update_ch = bm_updated(i,:);
    
    axes2 = axes('Parent',h,'Position',pv_bm);
    hold(axes2,'on');
%     ylim(axes2,[min(data)-.5 max(data)+.5]);
    xlim(axes2,[1 iter_num]);
    
    plot(1:iter_num, bm_update_ch, 'Color', color_bm, 'LineWidth', 1.5),hold on
    plot(1:iter_num, bm_real(i).*ones(1,iter_num), 'b', 'LineWidth', 2)
    title('b_m')
%     legend('estimated value', 'real value')
    
    hold on
    p = plot(axes2 , iter, bm_update_ch(iter) , 'o');
    p.MarkerSize = color_size;
    p.MarkerFaceColor = color_bm;
    p.MarkerEdgeColor = color_bm;
    
end

end

% plot mu
function [] = plot_mu(mu_estimate, mu_real, iter_num, pv_mu1, pv_mu2, pv_mu3, pv_mu4, h, iter, marker_size, color_real, color_imag)



mu_1_real = real(mu_estimate(1,:));
mu_1_imag = imag(mu_estimate(1,:));

mu_2_real = real(mu_estimate(2,:));
mu_2_imag = imag(mu_estimate(2,:));

%% 1
str_tit = 'real part of 1st element - \mu';
true_value = real(mu_real(1));
plot_mu_sub(mu_1_real, true_value, iter_num, h, pv_mu1, iter, str_tit, marker_size, color_real)

%% 2
str_tit = 'imaginary part of 1st element - \mu';
true_value = imag(mu_real(1));
plot_mu_sub(mu_1_imag, true_value, iter_num, h, pv_mu2, iter, str_tit, marker_size, color_imag)

%% 3
str_tit = 'real part of 2nd element - \mu';
true_value = real(mu_real(2));
plot_mu_sub(mu_2_real, true_value, iter_num, h, pv_mu3, iter, str_tit, marker_size, color_real)

%% 4
str_tit = 'imaginary part of 2nd element - \mu';
true_value = imag(mu_real(2));
plot_mu_sub(mu_2_imag, true_value, iter_num, h, pv_mu4, iter, str_tit, marker_size, color_imag)

end

function [] = plot_mu_sub(data, true_value, iter_num, h, positionVector, iter, str_tit, color_size, color_1)

axes1 = axes('Parent',h,'Position',positionVector);
hold(axes1,'on');
ylim(axes1,[min(data)-.5 max(data)+.5]);
xlim(axes1,[1 iter_num]);

plot(1:iter_num, data, 'LineWidth', 1.5, 'Color', color_1), 

hold on, plot(true_value*ones(1, length(data)), 'b')
title(str_tit)
xlabel('iter')
% legend('updated values', 'real value')

hold on
p = plot(axes1 , iter, data(iter) , 'o');
p.MarkerSize = color_size;
p.MarkerFaceColor = color_1;
p.MarkerEdgeColor = color_1;


end

% plot L
function [] = plot_L(L_estimate, L_real, pv_real_L, pv_imag_L, h, iter )

tit_real = sprintf('real(L_{esitmated} - L_{true}) - iter %d', iter);
L_real_diff = real(L_estimate(:,:,iter) - L_real);
plot_L_sub(L_real_diff, pv_real_L, h, tit_real)


tit_real = sprintf('imag(L_{esitmated} - L_{true}) - iter %d', iter);
L_imag_diff = imag(L_estimate(:,:,iter) - L_real);
plot_L_sub(L_imag_diff, pv_imag_L, h, tit_real);

end

function [] = plot_L_sub(L_diff, pv_L, h, tit)

axes1 = axes('Parent',h, 'Position',pv_L);
hold(axes1,'on');

imagesc(L_diff), hold on

title(tit)

end

% plot sv^2
function [] = plot_sv(sv_estimate, sv_real, iter_num)

figure, 
plot(1:iter_num, sv_estimate), hold on
plot(1:iter_num, sv_real*ones(1,iter_num), 'r', 'LineWidth',2)
title('sv^2 estimate')
legend('estimated sequence', 'real time')


end
