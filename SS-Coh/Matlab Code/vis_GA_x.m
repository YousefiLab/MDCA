function [] = vis_GA_x(Bayes_str, xk_sample, dim, iter_plot, strct_min_max, using_sim_data, t_min)

MU_filter = Bayes_str.MU_filter;
SIGMA_filter = Bayes_str.SIGMA_filter;
MU_smoother = Bayes_str.MU_smoother;
SIGMA_smoother = Bayes_str.SIGMA_smoother;


switch using_sim_data
    
    case 0
        fun_plt_cnfdnc_intvl_real(MU_smoother, SIGMA_smoother, 'Smoother Estimation', dim, iter_plot, strct_min_max, t_min)
        
    case 1
        %% plot seperately
        % fun_plt_cnfdnc_intvl(xk_sample, MU_filter, SIGMA_filter, 'Filter Estimation', dim, iter_plot)
        fun_plt_cnfdnc_intvl_sim(xk_sample, MU_smoother, SIGMA_smoother, 'Smoother Estimation', dim, iter_plot, strct_min_max, t_min)
        
        
        %% subplot of filter smoother
        % fun_subPlot_cnfdnc_intv(xk_sample, MU_filter, SIGMA_filter, MU_smoother, SIGMA_smoother, dim, iter_plot)
        
        %% plot summation
        if dim==2
            str_tit_sum = sprintf('summation of x_1 and x_2 --- %d Iter', iter_plot);
            % fun_plot_sum_2d(xk_sample,MU_filter, MU_smoother, str_tit_sum)
        end
        
        
        
end

end

%% nested functions 
function [] = fun_plt_cnfdnc_intvl_sim(x_sim, mu_estimation, sigma_estimation, tit, dim, iter_plot, strct_min_max, t_min)



m = size(mu_estimation);

cnfdnc_intvl = zeros(2*m(1), m(2));
% 1st row: up confidence interval of 1d
% 2nd row: down confidence interval of 1d
% 3rd row: up confidence interval of 2d
% 4th row: down confidence interval of 2d

for i=1: m(1) % number of dimensiona (2)
    
    for j=1: m(2) % number of points (127)
        mu_temp_dim = mu_estimation(i, j);
        sigma_temp_dim = squeeze(sigma_estimation(j, :, :));
        std_temp = sqrt(sigma_temp_dim(i, i));
        
        coef_ = 4;
        up_lim_temp = mu_temp_dim + coef_*std_temp;
        down_lim_temp = mu_temp_dim - coef_*std_temp;
        
        cnfdnc_intvl( (i-1)*dim + 1, j ) = up_lim_temp;
        cnfdnc_intvl( (i-1)*dim + 2, j ) = down_lim_temp;
    end
end

x_ind = 1: length(cnfdnc_intvl);

switch dim 
    case 1 
        min_1_d1 = strct_min_max.min_1_d1;
        max_1_d1 = strct_min_max.max_1_d1;
        dev_val = 0.2;
        
        cnfdnc_intvl_d1 = cnfdnc_intvl(1:2, :);
        
        % call plot function 
        funSub_plt_cnfdnt_sim(x_sim, x_ind, cnfdnc_intvl_d1, mu_estimation, tit, 1, iter_plot, min_1_d1-dev_val, max_1_d1+dev_val, t_min)
%         funSub_plt_cnfdnt_sim_igmSC(x_sim, x_ind, cnfdnc_intvl_d1, mu_estimation, tit, 1, iter_plot, min_1_d1-dev_val, max_1_d1+dev_val, t_min)
    case 2 
        min_1_d1 = strct_min_max.min_1_d1;
        max_1_d1 = strct_min_max.max_1_d1;
        min_1_d2 = strct_min_max.min_1_d2;
        max_1_d2 = strct_min_max.max_1_d2;
        dev_val = 0.1;
        
        
        cnfdnc_intvl_d1 = cnfdnc_intvl(1:2, :);
        cnfdnc_intvl_d2 = cnfdnc_intvl(3:4, :);
        
        % call plot function 
        funSub_plt_cnfdnt_sim(x_sim, x_ind, cnfdnc_intvl_d1, mu_estimation(1,:), x_sim(1,:), tit, 1, iter_plot, min_1_d1-dev_val, max_1_d1+dev_val, t_min)
        funSub_plt_cnfdnt_sim(x_sim, x_ind, cnfdnc_intvl_d2, mu_estimation(2,:), x_sim(2,:), tit, 2, iter_plot, min_1_d2-dev_val, max_1_d2+dev_val, t_min)
        
        
        figure, 
        hold all 
        plot(x_ind, x_sim(1,:), '--',  'LineWidth', 1, 'Color', [16, 50, 196]./255)
        plot(x_ind, x_sim(2,:), '--', 'LineWidth', 1, 'Color', [201, 24, 10]./255)
        plot(x_ind, mu_estimation(1,:), 'LineWidth', 2, 'Color', [16, 50, 196]./255)
        plot(x_ind, mu_estimation(2,:), 'LineWidth', 2, 'Color', [201, 24, 101]./255)
        legend('x_{d1}', 'x_{d2}', 'estimated /mu_1', 'estimated /mu_2')
        xlim([x_ind(1) x_ind(end)])
        
        str_tit = sprintf('%s -- %d iter', tit, iter_plot);
        title(str_tit)
        hold off 
end


end

function [] = fun_plt_cnfdnc_intvl_real(mu_estimation, sigma_estimation, tit, dim, iter_plot, strct_min_max, t_min)



m = size(mu_estimation);

cnfdnc_intvl = zeros(2*m(1), m(2));
% 1st row: up confidence interval of 1d
% 2nd row: down confidence interval of 1d
% 3rd row: up confidence interval of 2d
% 4th row: down confidence interval of 2d

for i=1: m(1) % number of dimensiona (2)
    
    for j=1: m(2) % number of points (127)
        mu_temp_dim = mu_estimation(i, j);
        sigma_temp_dim = squeeze(sigma_estimation(j, :, :));
        std_temp = sqrt(sigma_temp_dim(i, i));
        coef_ = 4; 
        up_lim_temp = mu_temp_dim + coef_*std_temp;
        down_lim_temp = mu_temp_dim - coef_*std_temp;
        
        cnfdnc_intvl( (i-1)*dim + 1, j ) = up_lim_temp+ .5;
        cnfdnc_intvl( (i-1)*dim + 2, j ) = down_lim_temp- .5;
    end
end

x_ind = 1: length(cnfdnc_intvl);

switch dim 
    case 1 
        min_1_d1 = strct_min_max.min_1_d1;
        max_1_d1 = strct_min_max.max_1_d1;
        dev_val = 0.2;
        
        cnfdnc_intvl_d1 = cnfdnc_intvl(1:2, :);
        
        % call plot function 
% % %         funSub_plt_cnfdnt(x_ind, cnfdnc_intvl_d1, mu_estimation, x_sim, tit, 1, iter_plot, min_1_d1-dev_val, max_1_d1+dev_val)
        funSub_plt_cnfdnt_real(x_ind, cnfdnc_intvl_d1, mu_estimation, tit, 1, iter_plot, min_1_d1-dev_val, max_1_d1+dev_val, t_min)
    case 2 
        min_1_d1 = strct_min_max.min_1_d1;
        max_1_d1 = strct_min_max.max_1_d1;
        min_1_d2 = strct_min_max.min_1_d2;
        max_1_d2 = strct_min_max.max_1_d2;
        dev_val = 0.1;
        
        
        cnfdnc_intvl_d1 = cnfdnc_intvl(1:2, :);
        cnfdnc_intvl_d2 = cnfdnc_intvl(3:4, :);
        
        % call plot function 
        funSub_plt_cnfdnt_real(x_ind, cnfdnc_intvl_d1, mu_estimation(1,:), x_sim(1,:), tit, 1, iter_plot, min_1_d1-dev_val, max_1_d1+dev_val, t_min)
        funSub_plt_cnfdnt_real(x_ind, cnfdnc_intvl_d2, mu_estimation(2,:), x_sim(2,:), tit, 2, iter_plot, min_1_d2-dev_val, max_1_d2+dev_val, t_min)
        
        
        figure, 
        hold all 
        plot(x_ind, x_sim(1,:), '--',  'LineWidth', 1, 'Color', [16, 50, 196]./255)
        plot(x_ind, x_sim(2,:), '--', 'LineWidth', 1, 'Color', [201, 24, 10]./255)
        plot(x_ind, mu_estimation(1,:), 'LineWidth', 2, 'Color', [16, 50, 196]./255)
        plot(x_ind, mu_estimation(2,:), 'LineWidth', 2, 'Color', [201, 24, 101]./255)
        legend('x_{d1}', 'x_{d2}', 'estimated /mu_1', 'estimated /mu_2')
        xlim([x_ind(1) x_ind(end)])
        
        str_tit = sprintf('%s -- %d iter', tit, iter_plot);
        title(str_tit)
        hold off 
end


end

function [] = funSub_plt_cnfdnt_real(x_ind, cnfdnc_intvl_data, mu_estimation, tit, dim_num, iter_plot, y_min, y_max, t_min)
fnt_size = 20; 
lgd_size = 12;
num_stp = 3;

step_x = t_min(end)./length(x_ind);
x_ind_t = 0: step_x: t_min(end)-step_x;

patch_color = [161, 204, 247; 161, 204, 247]./255;
mu_color = [17, 128, 240; 17, 128, 240]./255; 
x_color = [230, 23, 81 ;18, 17, 18]./255;

switch tit
    case 'Filter Estimation'
        patch_color = patch_color(1, :);
        mu_color = mu_color(1, :);
        
    case 'Smoother Estimation'
        patch_color = patch_color(2, :);
        mu_color = mu_color(2, :);
end

switch dim_num
    case 1
        x_color = x_color(1, :);
    case 2 
        x_color = x_color(2, :);
end



% h = figure('units','normalized','outerposition',[0 0 1 1]);
h = figure;

hold all
box on;
% cnfdnc_intvl_data = 5*cnfdnc_intvl_data;
% mu_estimation = 5*mu_estimation;
patch([x_ind_t fliplr(x_ind_t)], [cnfdnc_intvl_data(1,:) fliplr(cnfdnc_intvl_data(2,:))], patch_color, 'LineStyle','none')

plot(x_ind_t, mu_estimation(1,:), 'Color', mu_color, 'LineWidth', 1);


xlim([0 x_ind_t(end)])
str_tit = sprintf(' dim: %d -- Iter: %d', dim_num, iter_plot);
% title(str_tit)
xlabel('Time (min)')
ylabel({'Unconsciousness'; 'State (x)'})
% ylabel('Unconsciousness State ({\boldmath$x$})', 'Interpreter','latex', 'FontSize', fnt_size)
set(gca,'ytick',[])
set(gca,'yticklabel',[])
% ylim([y_min y_max])

ax = gca;

x_num_step = num_stp;
y_num_step = num_stp;

% xtick
x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested function 
ax.XTick = x_tick_values;

% ytick
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested function 
ax.YTick = y_tick_values;

ax.FontSize = fnt_size;

lgd = legend('Confidence Interval', 'Estimated $x_{k}$', 'interpreter','latex');
lgd.FontSize = lgd_size;
hold off

ax.FontSize = 28;

str_save = sprintf('/Results/x_smoother_real_%dd_iter%d.png', dim_num, iter_plot);
saveas(h,[pwd str_save]);

end

function [] = funSub_plt_cnfdnt_sim(x_sim, x_ind, cnfdnc_intvl_data, mu_estimation, tit, dim_num, iter_plot, y_min, y_max, t_min)
fnt_size = 20; 
lgd_size = 12;
num_stp = 3;

patch_color = [ 193, 198, 219; 184, 222, 192]./255;
mu_color = [ 5, 8, 176; 26, 176, 56]./255; 
x_color = [230, 23, 81 ;18, 17, 18]./255;

switch tit
    case 'Filter Estimation'
        patch_color = patch_color(1, :);
        mu_color = mu_color(1, :);
        
    case 'Smoother Estimation'
        patch_color = patch_color(2, :);
        mu_color = mu_color(2, :);
end

switch dim_num
    case 1
        x_color = x_color(1, :);
    case 2 
        x_color = x_color(2, :);
end

MSE = mean((x_sim - mu_estimation).^2);


% h = figure('units','normalized','outerposition',[0 0 1 1]);
h = figure;
hold all
patch([x_ind fliplr(x_ind)], [cnfdnc_intvl_data(1,:) fliplr(cnfdnc_intvl_data(2,:))], patch_color, 'LineStyle','none')
plot(mu_estimation(1,:), '--', 'Color', mu_color, 'LineWidth', 2);
plot(x_sim(1, :), 'Color', x_color, 'LineWidth', 2);

% xlim([0 t_min(end)])
xlim([1 x_ind(end)])
str_tit = sprintf('MSE: %0.4f -- %dd -- Iter: %d', MSE, dim_num, iter_plot);
title(str_tit)
xlabel('Trial')

ylim([y_min y_max])

ax = gca;

x_num_step = num_stp;
y_num_step = num_stp;

% xtick
x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested function 
ax.XTick = x_tick_values;

% ytick
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested function 
ax.YTick = y_tick_values;

ax.FontSize = fnt_size;

lgd = legend('Confidence Interval', 'Estimated $x_{k}$', 'Simulated  $x_{k}$', 'interpreter','latex');
lgd.FontSize = lgd_size;


% ylabel('Consciousness State', 'FontSize', fnt_size)
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])

ylabel('{\boldmath$x$}', 'Interpreter','latex', 'FontSize', fnt_size+10)

hold off

str_save = sprintf('/Results/x_smoother_sim_%dd_iter%d.png', dim_num, iter_plot);
saveas(h,[pwd str_save]);

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
cnt = find_order(step_tick); % call find_order
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

