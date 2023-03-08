function [] = vis_GA_W(W_sim, W_iter, dim, using_sim_data, use_bias)
mw = size(W_iter);
Iter = mw(1)-1;

% plot Euclidean Distance and Absolute Error TOT
vis_W_plot_1(W_sim, W_iter, Iter, dim, using_sim_data)

% plot updating values of W over time
vis_W_plot_2(W_sim, W_iter, Iter, using_sim_data)

% plot imagesc figures
vis_W_plot_3(W_sim, W_iter, Iter, using_sim_data)
end


%% nested
% plot Euclidean Distance and Absolute Error TOT
function [] = vis_W_plot_1(W_sim, W_iter, Iter, dim, using_sim_data)
fnt_size = 45;
tick_num_step = 3;
line_width = 5;

W_prev = squeeze(W_iter(1:end-1, :, :));
W_update = squeeze(W_iter(2:end, :, :));


%% plot Absolute Error W iteration to the preivous
W_diff_prev_cur = W_diff_update(W_prev, W_update, Iter);

h = figure('units','normalized','outerposition',[0 0 1 1]);
plot(W_diff_prev_cur, 'LineWidth', line_width)
xlabel('iteratoin'), title('Absolute Error between $W_{current}$ and $W_{prev}$', 'interpreter', 'latex')
xlim([1 Iter])
ax = gca;

if Iter > 20
    % xtick
    x_num_step = tick_num_step;
    x_tick_values = fun_finding_tick(ax, x_num_step, 'x');
    ax.XTick = x_tick_values;
    
    % ytick
    y_num_step = tick_num_step;
    y_tick_values = fun_finding_tick(ax, y_num_step, 'y');
    ax.YTick = y_tick_values;
end

ax.FontSize = fnt_size;
str_save = sprintf('/Results/W_%dd_absError_iter%d.png', dim, Iter);
saveas(h,[pwd str_save]);

%% plot Euclidean Distance
if using_sim_data==1
    euclidean_dist = fun_euclidean_dist(W_sim, W_prev, Iter); % nested fun
    
    mm = size(W_sim);
    
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(euclidean_dist, 'Linewidth', line_width)
    xlabel('iteratoin'), title('Euclidean Distance $W_{simulation}$ and $W_{training}$', 'interpreter', 'latex')
    xlim([1 Iter])
    ax = gca;
    
    if Iter > 20
        % xtick
        x_num_step = tick_num_step;
        x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
        ax.XTick = x_tick_values;
        
        % ytick
        y_num_step = tick_num_step;
        y_tick_values = fun_finding_tick(ax, y_num_step, 'y');% nested fun
        ax.YTick = y_tick_values;
    end
    
    ax.FontSize = fnt_size;
    str_save = sprintf('/Results/W_%dd_EuclideanDist_iter%d.png', dim, Iter);
    saveas(h,[pwd str_save]);
end

end

% plot updating values of W over time
function [] = vis_W_plot_2(W_sim, W_all_itrtn, Iter, using_sim_data)

sim_color = [0 0 0]./255;
training_color = [19, 208, 212; 230, 44, 161; 138, 45, 101]./255;
line_width_sim = 3;
line_width_training = 2;
fnt_size = 25;
ylim_bias_val = 1;
tick_num_step = 3;

m = size(W_sim);


for col_dim=1: m(2) % dimension+1
    
    for row_ch=1: m(1) % num of channel
        w_sim_vec = W_sim(row_ch, col_dim)*ones(1, Iter+1);
        h = figure;
        if using_sim_data==1
            plot(w_sim_vec, 'LineWidth', line_width_sim, 'Color', sim_color)
        end
        hold on
        w_training_vec = fun_get_updated_vec(W_all_itrtn, row_ch, col_dim);
        plot(w_training_vec, '--', 'LineWidth', line_width_training, 'Color', training_color(col_dim, :));
        str_tit = sprintf('col: %d - row(channel): %d', col_dim, row_ch);
        title(str_tit)
        if using_sim_data==1
            legend('W simulation', 'W training', 'Interpreter', 'latex')
        end
        
        max_val = max(max(w_sim_vec), max(w_training_vec));
        min_val = min(min(w_sim_vec), min(w_training_vec));
        
        ylim([min_val-ylim_bias_val max_val+ylim_bias_val])
        xlim([1 Iter+1])
        
        ax = gca;
        if Iter > 20
            % xtick
            x_num_step = tick_num_step;
            x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
            ax.XTick = x_tick_values;
            
            % ytick
            y_num_step = tick_num_step;
            y_tick_values = fun_finding_tick(ax, y_num_step, 'y');% nested fun
            ax.YTick = y_tick_values;
        end
        
        ax.FontSize = fnt_size;
        str_save = sprintf('/Results/W_%dd_col%d_row%d_iter%d.png', m(2)-1, col_dim, row_ch, Iter);
        saveas(h,[pwd str_save]);
        
        
    end
end



end

% plot imagesc figures
function [] = vis_W_plot_3(W_sim, W_all_itrtn, Iter, using_sim_data)
fnt_size = 20;

mm =size(W_sim);
y_tick_vec = 1: mm(1);
x_tick_vec = 1: mm(2);

W_init = squeeze(W_all_itrtn(1, :, :));
W_final = squeeze(W_all_itrtn(Iter+1, :, :));


%% plot W
if using_sim_data==0
h = figure;
imagesc(W_init.')
xticks(y_tick_vec)
yticks(x_tick_vec)
yticklabels({'slope', 'bias'})

ylabel('Channel')
ax = gca;

title('W init')
colormap parula
cbh = colorbar;

min_val = ceil(min(W_init(:)));
max_val = floor(max(W_init(:)));

min_tot = min(abs(min_val), max_val);
if min_tot==0
    min_tot = max_val;
end

c_tick_vec = [-min_tot 0 min_tot];
cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(c_tick_vec);

ax.FontSize = fnt_size;
ylabel('{\boldmath$W$}','Interpreter','latex', 'FontSize', fnt_size+ 10)



str_save = sprintf('/Results/W_init_%dd_W_iter%d.png', mm(2)-1, Iter);
saveas(h,[pwd str_save]);

h=figure;
imagesc(W_final.')
xticks(y_tick_vec)
yticks(x_tick_vec)
yticklabels({'slope', 'bias'})

ylabel('Channel')
ax = gca;
ax.FontSize = fnt_size;
title('W final')

colormap parula
cbh = colorbar;

min_val = ceil(min(W_final(:)));
max_val = floor(max(W_final(:)));

min_tot = min(abs(min_val), max_val);
if min_tot==0
    min_tot = max_val;
end
c_tick_vec = [-min_tot 0 min_tot];
cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(c_tick_vec);

max_tot = max(abs(min(W_final(:))), max(W_final(:)));
caxis([-max_tot max_tot])

ax.FontSize = fnt_size;
ylabel('{\boldmath$W$}','Interpreter','latex', 'FontSize', fnt_size+ 10)

str_save = sprintf('/Results/W_fin_%dd_W_iter%d.png', mm(2)-1, Iter);
saveas(h,[pwd str_save]);

end


%% plot of W-Wsim
if using_sim_data==1
    W_init_sim = W_init - W_sim;
    W_final_sim = W_final - W_sim;
    
    min_val = min(min(W_final_sim(:)), min(W_init_sim(:)));
    max_val = max(max(W_final_sim(:)), max(W_init_sim(:)));
    
    h = figure;
    imagesc(W_init_sim.')
    xticks(y_tick_vec)
    yticks(x_tick_vec)
    yticklabels({'slope', 'bias'})
    
    xlabel('Channel')
    ax = gca;
    %     title('init')
    title('W init  - W sim')
    
    
    colormap(parula) ; %Create Colormap
    caxis([min_val max_val]);
    cbh = colorbar ; %Create Colorbar
    c_tick_vec = [ceil(min_val) 0 floor(max_val)];
    cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
    cbh.TickLabels = num2cell(c_tick_vec);
    ax.FontSize = fnt_size;
    
    ylabel('{\boldmath${dW}$}','Interpreter','latex', 'FontSize', fnt_size+ 10)
    
    str_save = sprintf('/Results/Winit_%dd_W_Wsim_iter%d.png', mm(2)-1, Iter);
    saveas(h,[pwd str_save]);
    
    %% final
    
    % setup values
    left_pos = .25;
    bot_pos = .25;
    width_val = .7;
    height_val = .33;
    
    
    bot_diff = .37;
    
    h=figure;
    % dW
    ax = subplot(212);
    ax.Position = [left_pos (bot_pos+ 0*bot_diff) width_val height_val];
    imagesc(W_final_sim.')
    xticks(y_tick_vec)
    yticks(x_tick_vec)
    yticklabels({'slope', 'bias'})
    
    
    xlabel('Channel')%     title('final')
%     title('W final  - W sim')
    
        
    colormap(parula) ; %Create Colormap
    caxis([min_val max_val]);
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
    cbh.TickLabels = num2cell(c_tick_vec);
    
    ax.FontSize = fnt_size;
    
    ylabel('{\boldmath${d}$}{\bf{W}}','Interpreter','latex', 'FontSize', fnt_size+ 10)
    
    
    % w final 
    ax = subplot(211);
    ax.Position = [left_pos (bot_pos+ 1*bot_diff) width_val height_val];
    
    imagesc(W_final.')
    xticks(y_tick_vec)
    yticks(x_tick_vec)
    yticklabels({'slope', 'bias'})
    
    ylabel('Channel')
    ax = gca;
    ax.FontSize = fnt_size;
%     title('W final')
    
    colormap parula
    cbh = colorbar;
    
    min_val = ceil(min(W_final(:)));
    max_val = floor(max(W_final(:)));
    
    min_tot = min(abs(min_val), max_val);
    if min_tot==0
        min_tot = max_val;
    end
    c_tick_vec = [-min_tot 0 min_tot];
    cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
    cbh.TickLabels = num2cell(c_tick_vec);
    
    set(ax,'xtick',[])
    set(ax,'xticklabel',[])
    
    max_tot = max(abs(min(W_final(:))), max(W_final(:)));
    caxis([-max_tot max_tot])
    
    ax.FontSize = fnt_size;
    ylabel('{\bf{W}}','Interpreter','latex', 'FontSize', fnt_size+ 10)
%     
%     str_save = sprintf('/Results/W_fin_%dd_W_iter%d.png', mm(2)-1, Iter);
%     saveas(h,[pwd str_save]);
    
%     h = figure;
    
    
    str_save = sprintf('/Results/Wfinal_%dd_W_Wsim_iter%d.png', mm(2)-1, Iter);
    saveas(h,[pwd str_save]);
    
end

end

% fun for Euclidean Distance
function euclidean_dist = fun_euclidean_dist(W_sim, W_prev, Iter)

euclidean_dist = zeros(1, Iter);
for iter=1: Iter
    w_prev_temp = squeeze(W_prev(iter, :, :));
    W_diff_prev = W_sim - w_prev_temp;
    
    W_diff_prev_vec = W_diff_prev(:);
    
    vec_diff = W_diff_prev_vec';
    D = sqrt( vec_diff* vec_diff');
    euclidean_dist(1, iter) = D;
end

end

% finding w training value for each element of W
function w_train_vec = fun_get_updated_vec(W_all_itrtn, ch_num, col_dim)
m = size(W_all_itrtn);

w_train_vec = zeros(1, m(1));
for iter=1: m(1)
    w_train_vec(1, iter) = squeeze(W_all_itrtn(iter,  ch_num, col_dim));
end

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


% find absolute differnce between 2 step of W
function [W_diff_prev_cur] = W_diff_update(W_prev, W_update, Iter)

W_diff_prev_cur = zeros(1, Iter);
for i=1:Iter
    
    W_diff_prev_temp = squeeze(W_prev(i, :, :));
    W_diff_prev_vec = W_diff_prev_temp(:);
    W_diff_update_temp = squeeze(W_update(i, :, :));
    W_diff_update_vec = W_diff_update_temp(:);
    
    
    W_diff_prev_cur(1, i) = sum(abs(W_diff_prev_vec-W_diff_update_vec));
end

end