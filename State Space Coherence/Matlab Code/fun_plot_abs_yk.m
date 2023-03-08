function [] = fun_plot_abs_yk(Yk, using_sim, dim)

point_each_plt = 20;
tick_num_step = 3;
fnt_size =20; 

my = size(Yk);
x = (1: my(1));

abs_yk = abs(Yk);

% figure, imagesc(abs_yk);

h=figure;
imagesc(abs_yk.');
set(gca,'YDir','normal')
colormap(parula)


xlabel('Trial')
set(gca,'YTickLabel',[]);


ax = gca;
x_num_step = tick_num_step;
x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
ax.XTick = x_tick_values;
colorbar
ax = gca;
ax.FontSize = fnt_size;

max_val = max(abs_yk(:));
c_tick_vec = [0 floor(max_val)./2 floor(max_val)];

colormap(parula) ; %Create Colormap
caxis([0 max_val]);
cbh = colorbar ; %Create Colorbar
cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(c_tick_vec);

ylabel('{\bf{Y}}','Interpreter','latex', 'FontSize', fnt_size+5)

if using_sim ==1
    str_save = sprintf('/Results/y_sim_%dd.png', dim);
elseif using_sim==0
    str_save = sprintf('/Results/y_real_%dd.png', dim);
end


saveas(h,[pwd str_save]);



aa = 0;
while aa==1
    aa=0;
    
    max_val = ceil(max(abs_yk(:)));
    
    step_y = max_val/point_each_plt;
    y_step = 0: step_y: max_val-step_y;
    
    mat_imgSC = zeros(my(1), point_each_plt*my(2));
    for i=1: my(2)
        temp_imgSC = zeros(my(1), point_each_plt);
        for j=1: my(1)
            
            abs_yk_temp = abs_yk(j, i);
            
            [val,idx] = min(abs(y_step-abs_yk_temp));
            %     minVal=y_step(idx);
            temp_imgSC(j, idx) = 1;
        end
        
        low_x = (i-1)*point_each_plt+1;
        up_x = i*point_each_plt;
        
        mat_imgSC(:, low_x: up_x) = temp_imgSC;
        
        
    end
    
    ax = figure;
    imagesc(mat_imgSC.')
    set(gca,'YDir','normal')
    xlabel('Trial')
    set(gca,'YTickLabel',[]);
    ylabel('Y')
    
    ax = gca;
    x_num_step = tick_num_step;
    x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
    ax.XTick = x_tick_values;
    
end

bb = 0;
while bb==1
    
    y_lim = [-1 max_val+1];
    
    tick_num_step = 3;
    y_ticks = [0 4];
    %% setup values
    left_pos = .1;
    bot_pos = .15;
    width_val = .7;
    height_val = .13;
    
    
    bot_diff = .14;
    
    line_width = .2;
    fnt_size= 10;
    
    %% plots
    h =figure;
    ax = subplot(5,1,5);
    plot(x, abs_yk(:, 1), 'LineWidth', line_width)
    ax.Position = [left_pos (bot_pos+ 0*bot_diff) width_val height_val];
    ylabel('ch 1')
    xlabel('Trial')
    yticks(y_ticks)
    ylim(y_lim)
    ax.FontSize = fnt_size;
    % xtick
    x_num_step = tick_num_step;
    x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
    ax.XTick = x_tick_values;
    
    ax.XAxis.FontSize = fnt_size + 8;
    ax.XLabel.FontSize = fnt_size + 8;
    
    
    ax = subplot(5,1,4);
    plot(x, abs_yk(:, 2), 'LineWidth', line_width)
    ax.Position = [left_pos (bot_pos+ 1*bot_diff) width_val height_val];
    ylabel('ch 2')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    yticks(y_ticks)
    ylim(y_lim)
    ax.FontSize = fnt_size;
    
    ax = subplot(5,1,3);
    plot(x, abs_yk(:, 3), 'LineWidth', line_width)
    ax.Position = [left_pos (bot_pos+ 2*bot_diff) width_val height_val];
    ylabel('ch 3')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    yticks(y_ticks)
    ylim(y_lim)
    ax.FontSize = fnt_size;
    
    ax = subplot(5,1,2);
    plot(x, abs_yk(:, 4), 'LineWidth', line_width)
    ax.Position = [left_pos (bot_pos+ 3*bot_diff) width_val height_val];
    ylabel('ch 4')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    yticks(y_ticks)
    ylim(y_lim)
    ax.FontSize = fnt_size;
    
    ax = subplot(5,1,1);
    plot(x, abs_yk(:, 5), 'LineWidth', line_width)
    ax.Position = [left_pos (bot_pos+ 4*bot_diff) width_val height_val];
    ylabel('ch 5')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    yticks(y_ticks)
    ylim(y_lim)
    ax.FontSize = fnt_size;
    
    
    
    % [left bottom width height]
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
