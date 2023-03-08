function [] = vis_GA_mu(mu_sim, MU_iter, dim, using_sim_data)
fnt_size = 45;
tick_num_step = 3;
line_width = 5;

mu_sim = mu_sim.';

mw = size(MU_iter);
Iter = mw(1) -1;

mu_prev = (MU_iter(1:end-1, :));
mu_update = (MU_iter(2:end, :));

[abs_mu_diff_prev_cur, abs_mu_diff_prev_sim] = mu_diff_update(mu_prev, mu_update, mu_sim, Iter);

plt_like_w = 0
if plt_like_w == 1
% new part 
    real_mu_sim = real(mu_sim);
    imag_mu_sim = imag(mu_sim);
    MU_SIM = zeros(length(real_mu_sim), 2);
    MU_SIM(:, 1) = real_mu_sim';
    MU_SIM(:, 2) = imag_mu_sim';

    real_mu_iter = real(MU_iter);
    image_mu_iter = imag(MU_iter);

    mm = size(MU_iter);
    MU_ALL_ITER = zeros(mm(1), mm(2), 2);
    for ii=1: mm(1)
        MU_ALL_ITER(ii, :, :) = [real_mu_iter(ii, :)' image_mu_iter(ii, :).'];
    end

    vis_W_plot_3(MU_SIM, MU_ALL_ITER, Iter, using_sim_data)
end


%% plot Absolute Error of mu - current iteration to the preivous
h = figure('units','normalized','outerposition',[0 0 1 1]);
plot(abs_mu_diff_prev_cur, 'LineWidth', line_width)
xlabel('iteratoin'), title('Absolute Error between $\mu_{current}$ and $\mu_{prev}$', 'interpreter', 'latex')
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
str_save = sprintf('/Results/mu_%dd_absError_prev_curr_iter%d.png', dim, Iter);
saveas(h,[pwd str_save]);


%% plot Absolute Error of mu - current iteration to the preivous
if using_sim_data == 1
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(abs_mu_diff_prev_sim, 'LineWidth', line_width)
    xlabel('iteratoin'), title('Absolute Error between $\mu_{current}$ and $\mu_{sim}$', 'interpreter', 'latex')
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
    str_save = sprintf('/Results/mu_%dd_absError_prev_sim_iter%d.png', dim, Iter);
    saveas(h,[pwd str_save]);
    
end
end

%% nested functions
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
ylabel('{\boldmath$\mu$}','Interpreter','latex', 'FontSize', fnt_size+ 10)



str_save = sprintf('/Results/MU_init_%dd_MU_iter%d.png', mm(2)-1, Iter);
saveas(h,[pwd str_save]);

h=figure;
imagesc(W_final.')
xticks(y_tick_vec)
yticks(x_tick_vec)
yticklabels({'Real', 'Imag'})

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
ylabel('{\boldmath$\mu$}','Interpreter','latex', 'FontSize', fnt_size+ 10)

str_save = sprintf('/Results/MU_fin_%dd_MU_iter%d.png', mm(2)-1, Iter);
saveas(h,[pwd str_save]);

end


%% plot of W-Wsim
if using_sim_data==1
    W_init_sim = W_init - W_sim;
    W_final_sim = W_final - W_sim;
    
    min_val = min(min(W_final_sim(:)), min(W_init_sim(:)));
    max_val = max(max(W_final_sim(:)), max(W_init_sim(:)));
    
    show_single_mu=0;
    if show_single_mu ==1
        h = figure;
        imagesc(W_init_sim.')
        xticks(y_tick_vec)
        yticks(x_tick_vec)
        yticklabels({'Real', 'Imag'})
        
        xlabel('Channel')
        ax = gca;
        %     title('init')
        title('MU init  - MU sim')
        
        
        colormap(parula) ; %Create Colormap
        caxis([min_val max_val]);
        cbh = colorbar ; %Create Colorbar
        hAxes = gca;
        colourRange = caxis( hAxes );
        cbh.Ticks = linspace( colourRange(1), colourRange(end), 3 );
        
        %     min_tick = ceil(min_val);
        %     max_tick = floor(max_val);
        %     c_tick_vec = [min_tick mean([min_tick max_tick]) max_tick];
        %     cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
        %     cbh.TickLabels = num2cell(c_tick_vec);
        ax.FontSize = fnt_size;
        
        ylabel('{\boldmath${d\mu}$}','Interpreter','latex', 'FontSize', fnt_size+ 10)
        
        str_save = sprintf('/Results/MUinit_%dd_MU_MUsim_iter%d.png', mm(2)-1, Iter);
        saveas(h,[pwd str_save]);
    end
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
    yticklabels({'Real', 'Imag'})
    
    
    xlabel('Channel')%     title('final')
%     title('W final  - W sim')
    
        
    colormap(parula) ; %Create Colormap
%     caxis([min_val max_val]);
    cbh = colorbar ; %Create Colorbar
    hAxes = gca;
    colourRange = caxis( hAxes );
    cbh.Ticks = linspace( colourRange(1), colourRange(end), 3 );
%     cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
%     cbh.TickLabels = num2cell(c_tick_vec);
    
    ax.FontSize = fnt_size;
    
    ylabel('\boldmath${d\mu}$','Interpreter','latex', 'FontSize', fnt_size+ 10)
    
    
    % w final 
    ax = subplot(211);
    ax.Position = [left_pos (bot_pos+ 1*bot_diff) width_val height_val];
    
    imagesc(W_final.')
    xticks(y_tick_vec)
    yticks(x_tick_vec)
    yticklabels({'Real', 'Imag'})
    
    ylabel('Channel')
    ax = gca;
    ax.FontSize = fnt_size;
%     title('W final')
    
    colormap parula
    cbh = colorbar;
    
    hAxes = gca;
    colourRange = caxis( hAxes );
    cbh.Ticks = linspace( colourRange(1), colourRange(end), 3 );
    
    
%     min_val = ceil(min(W_final(:)));
%     max_val = floor(max(W_final(:)));
%     
%     min_tot = min(abs(min_val), max_val);
%     if min_tot==0
%         min_tot = max_val;
%     end
%     c_tick_vec = [-min_tot 0 min_tot];
%     cbh.Ticks =  c_tick_vec; %Create 8 ticks from zero to 1
%     cbh.TickLabels = num2cell(c_tick_vec);
    
    set(ax,'xtick',[])
    set(ax,'xticklabel',[])
    
    max_tot = max(abs(min(W_final(:))), max(W_final(:)));
    caxis([-max_tot max_tot])
    
    ax.FontSize = fnt_size;
    ylabel('\boldmath${\mu}$','Interpreter','latex', 'FontSize', fnt_size+ 10)
%     
%     str_save = sprintf('/Results/W_fin_%dd_W_iter%d.png', mm(2)-1, Iter);
%     saveas(h,[pwd str_save]);
    
%     h = figure;
    
    
    str_save = sprintf('/Results/MUfinal_%dd_MU_MUsim_iter%d.png', mm(2)-1, Iter);
    saveas(h,[pwd str_save]);
    
end

end




% find absolute differnce between 2 step of W
function [abs_mu_diff_prev_cur, abs_mu_diff_prev_sim] = mu_diff_update(mu_prev, mu_update, mu_sim, Iter)

abs_mu_diff_prev_cur = zeros(1, Iter);
abs_mu_diff_prev_sim = zeros(1, Iter);
for i=1:Iter
    
    mu_diff_prev_temp = (mu_prev(i, :));
    mu_diff_update_temp = (mu_update(i, :));
    
    vec_diff_curr_prev = mu_diff_prev_temp-mu_diff_update_temp;
    vec_diff_prev_sim = mu_diff_prev_temp-mu_sim;
    
    abs_mu_diff_prev_cur(1, i) = sum(abs(vec_diff_curr_prev));
    abs_mu_diff_prev_sim(1, i) = sum(abs(vec_diff_prev_sim));
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
