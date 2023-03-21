function [] = fun_plt_simulated_data(xk_sample, dim, strct_min_max)

d_val = .2;
fnt_size = 20;
tick_num_step = 3;
line_width = 2;

if dim==2
    
%     h = figure('units','normalized','outerposition',[0 0 1 1]);
    h = figure;
    x_color = [148, 28, 1; 18, 17, 18]./255;
    hold all
    plot(1:length(xk_sample), xk_sample(1,:), 'LineWidth', line_width, 'Color', x_color(1, :));
    plot(1:length(xk_sample), xk_sample(2,:), 'LineWidth', line_width, 'Color', x_color(2, :))
    xlim([1 length(xk_sample)]);
    title('simulated data')
    xlabel('trial')
    legend('1d', '2d')
    hold off
    
    ax = gca;
    % xtick
    x_num_step = tick_num_step;
    x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
    ax.XTick = x_tick_values;

    % ytick
    y_num_step = tick_num_step;
    y_tick_values = fun_finding_tick(ax, y_num_step, 'y');% nested fun
    ax.YTick = y_tick_values;
    
    ax.FontSize = fnt_size;
    
    str_save = sprintf('x_sim_%dd.png', dim);
    saveas(h, str_save)
    
elseif dim==1
    
    min_1_d1 = strct_min_max.min_1_d1;
    max_1_d1 = strct_min_max.max_1_d1;
    
    
%     h = figure('units','normalized','outerposition',[0 0 1 1]);
    h = figure;
    x_color = [230, 23, 81]./255;
    hold all
    plot(1:length(xk_sample), xk_sample(1,:), 'LineWidth', line_width, 'Color', x_color(1, :));
    xlim([1 length(xk_sample)]);
    ylim([min_1_d1-d_val max_1_d1+d_val])
    title('simulated data')
    xlabel('Trial')
    
    hold off
    
    ax = gca;
    % xtick
    x_num_step = tick_num_step;
    x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
    ax.XTick = x_tick_values;

    % ytick
    y_num_step = tick_num_step;
    y_tick_values = fun_finding_tick(ax, y_num_step, 'y');% nested fun
    ax.YTick = y_tick_values;
    
    ax.FontSize = fnt_size;
    ylabel('{\boldmath$x$}','Interpreter','latex', 'FontSize', fnt_size+ 10)
    
    str_save = sprintf('/Results/x_sim_%dd.png', dim);
    saveas(h,[pwd str_save]);

end

end

%% nested functions 
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
