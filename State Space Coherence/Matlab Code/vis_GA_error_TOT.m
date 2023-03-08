function [] = vis_GA_error_TOT(abs_error_TOT, dim, Iter)
fnt_size = 45;
fnt_ttl = 35;
tick_num_step = 3;
line_width = 5; 

Iter = length(abs_error_TOT);

h = figure('units','normalized','outerposition',[0 0 1 1]);
plot(abs_error_TOT, 'LineWidth', line_width)
xlabel('iteratoin'), 
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
title('Absolute Error between $param_{current}$ and $param_{prev}$', 'interpreter', 'latex', 'FontSize', fnt_ttl)


str_save = sprintf('/Results/absErrorTot_%dd_iter%d.png', dim, Iter);
saveas(h, [pwd str_save]);

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
