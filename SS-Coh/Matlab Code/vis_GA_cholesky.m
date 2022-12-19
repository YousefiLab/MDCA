function [] = vis_GA_cholesky(xk_sample, Yk_sample, cnfg_sim, W_iter, MU_iter, dim, bayes_first, bayes_last, Iter)

    
   w_type = 'init';
    x_tick_init = [-2, 0, 2];
    x_lim_init = [-4 4]; 
    y_tick_init = [-2, 0, 2];
    y_lim_init = [-4 4]; 
    W_init = squeeze(W_iter(1, :, :));
    x_inint = bayes_first.MU_smoother;
    mu_init = MU_iter(1, :).';
    % mu_init_real = cnfg_sim.mu;
    Yk_sample_init = Yk_sample;
    Y_mapped_chol_init = fun_calc_cholesky(x_inint, Yk_sample_init, W_init, cnfg_sim.L, mu_init, dim);
    vis_cholesky(Y_mapped_chol_init, w_type,  x_tick_init, x_lim_init, y_tick_init, y_lim_init)



    w_type = 'sim';
    x_tick_sim = [-2, 0, 2];
    x_lim_sim = [-4 4]; 
    y_tick_sim = [-2, 0, 2];
    y_lim_sim = [-4 4]; 
    Y_mapped_chol_sim = fun_calc_cholesky(xk_sample, Yk_sample, cnfg_sim.W, cnfg_sim.L, cnfg_sim.mu, dim);
    vis_cholesky(Y_mapped_chol_sim, w_type, x_tick_sim, x_lim_sim, y_tick_sim, y_lim_sim)
%     
    w_type = 'init';
    x_tick_init = [-2, 0, 2];
    x_lim_init = [-4 4]; 
    y_tick_init = [-2, 0, 2];
    y_lim_init = [-4 4]; 
    W_init = squeeze(W_iter(1, :, :));
    x_inint = bayes_first.MU_smoother;
    mu_init = MU_iter(1, :).';
    % mu_init_real = cnfg_sim.mu;
    Yk_sample_init = Yk_sample;
    Y_mapped_chol_init = fun_calc_cholesky(x_inint, Yk_sample_init, W_init, cnfg_sim.L, mu_init, dim);
    vis_cholesky(Y_mapped_chol_init, w_type,  x_tick_init, x_lim_init, y_tick_init, y_lim_init)
    
    w_type = 'final';
    x_tick_fin = [-2, 0, 2];
    x_lim_fin = [-4 4]; 
    y_tick_fin = [-2, 0, 2];
    y_lim_fin = [-4 4]; 
    W_final = squeeze(W_iter(Iter, :, :));
    x_final = bayes_last.MU_smoother;
    mu_final = MU_iter(Iter, :).';
    Yk_sample_final = Yk_sample;
    Y_mapped_chol_final = fun_calc_cholesky(x_final, Yk_sample_final, W_final, cnfg_sim.L, mu_final, dim);
    vis_cholesky(Y_mapped_chol_final, w_type,  x_tick_fin, x_lim_fin, y_tick_fin, y_lim_fin)


end

%% nested functions 
% calc cholesky decomposition 
function Y_mapped_chol = fun_calc_cholesky(x, Yk, W, L, mu, dim)

mm = size(W);
M = length(x);

my = size(Yk);
Y_mapped_chol = zeros(my(2), my(1));



for i=1 : M
    x_k = x(:, i);
    

    D_eval = gen_sigma(x_k, W, dim); % nested fun
    D_eval = inv(D_eval);

    Chol_mat = sqrt(D_eval)*L';
    y_k = Yk(i, :);
    yk_remMu = (y_k-mu.');

    Y_mapped_chol(:, i) = Chol_mat*yk_remMu.';
    
end


end

% calc sgima based on W and x
function Sigma_W_x = gen_sigma(x_k, W, dim)
mm = size(W);
switch dim
    
    %1 d
    case 1
        Sigma_W_x = zeros(mm(1));
        for j=1: mm(1)
            w_m = W(j, 1);
            w_bias = W(j, 2);
            lambda_m = exp(x_k(1,1) * w_m(1,1)+ w_bias);
            Sigma_W_x(j, j) = lambda_m;
        end
        
        % 2d
    case 2
        Sigma_W_x = zeros(mm(1));
        for j=1: mm(1)
            w_m = W(j, 1:2);
            w_bias = W(j,3);
            lambda_m = exp(x_k(1,1) * w_m(1,1)  +   x_k(2, 1) * w_m(1, 2) + w_bias);
            Sigma_W_x(j, j) = lambda_m;
        end
end

end

% visualize cholesky 
function [] = vis_cholesky(Y_mapped_chol, w_type, x_tick, x_lim, y_tick, y_lim)
nbins = 125;
line_width = 4;
tick_num_step = 5;
marker_size = 4;
fnt_size =15;

Y_mapped_chol = real(Y_mapped_chol);

my = size(Y_mapped_chol);

switch w_type
    case 'sim'
        str_tit_1 = 'W Simulation -- ';
        str_save_1 = 'W_sim';
        color_hist = [83, 148, 100]./255;
        color_line = [5, 245, 68]./255;
        
    case 'init'
        str_tit_1 = 'W Initial -- ';
        str_save_1 = 'W_init';
        color_hist = [135, 53, 104]./255;
        color_line = [245, 2, 153]./255;
        
    case 'final'
        str_tit_1 = 'W Final -- ';
        str_save_1 = 'W_final';
        color_hist = [73, 152, 158]./255;
        color_line = [12, 224, 240]./255;
end

% plot 2d histograms
plt_2d = 1;
while plt_2d == 1
    plt_2d = 0;
for ch_num_1 =1: my(1)
    for ch_num_2 = ch_num_1+1: my(1)
        ch_1 = Y_mapped_chol(ch_num_1, :);
        ch_2 = Y_mapped_chol(ch_num_2, :);
        
        min_ = floor(min([ch_1, ch_2]));
        max_ = ceil(max([ch_1, ch_2]));
        
        max_t = max([abs(min_), max_]);
        x_lim_ = [-max_t, max_t];
        
        hh = figure('Renderer', 'painters', 'Position', [50 50 590 600]);
        sh = scatterhist(ch_1, ch_2,'Kernel','on', 'Color', color_hist);
        
        % label of scatter
        ax1 = sh(1);
        str_ch1 = sprintf('Ch %d', ch_num_1);
        str_ch2 = sprintf('Ch %d', ch_num_2);
        set(get(ax1,'XLabel'),'String', str_ch1)
        set(get(ax1,'YLabel'),'String', str_ch2)
        set(ax1,'XLim', x_lim_);
        set(ax1,'YLim', x_lim_);
%         set(ax1,'XTick', x_tick);
%         set(ax1,'YTick', y_tick);
        
        set(ax1,'FontSize', fnt_size);
        
%         ax.XTick = x_tick_values;
        
        % scatter settign
        x1=get(gca,'children');
        
        set(x1,'markerfacecolor', color_hist)
        set(x1,'markeredgecolor', color_hist)
        set(x1,'markersize', marker_size)
        
        % histogram setting
%         sh(2).Children(1).FaceColor = color_line;
        sh(2).Children(1).Color = 'none';
%         sh(2).Children(1).NumBins = nbins;
%         sh(3).Children(1).FaceColor = color_line;
        sh(3).Children(1).Color = 'none';
%         sh(3).Children(1).NumBins = nbins;
        
        str_tit = sprintf('%s - Joint Plot of  Ch #%d & Ch #%d', str_tit_1, ch_num_1, ch_num_2);
        suptitle(str_tit)
        
        str_save = sprintf('/Results/mutual_%s_ch%d_%d.png', str_save_1, ch_num_1, ch_num_2);
        saveas(hh,[pwd str_save]);
    end
end
end

% plot 1d histogram
plt_1d = 1;
while plt_1d == 1
    plt_1d = 0;
    for ch_num =1: my(1)
        
        y_ch = Y_mapped_chol(ch_num, :);
        
        pd = fitdist(y_ch.','Normal');
        
        min_val = min(y_ch);
        max_val = max(y_ch);
        
        max_abs = max(abs(min_val), abs(max_val));
        
        x_pdf = [-max_abs:0.1:max_abs];
        y_pdf = pdf(pd, x_pdf);
        
        hh = figure;
        h = histogram(y_ch, 'Normalization','pdf');
        h.NumBins = nbins;
        h.FaceColor = color_hist;
        h.EdgeColor = color_hist;
        
        l = line(x_pdf, y_pdf);
        l.LineWidth = line_width;
        l.Color = color_line;
        
        ax = gca;
        % xtick
        %     x_num_step = tick_num_step;
        %     x_tick_values = fun_finding_tick(ax, x_num_step, 'x');
        %     ax.XTick = x_tick_values;
        
        % ytick
        %     y_num_step = tick_num_step-2;
        %     y_tick_values = fun_finding_tick(ax, y_num_step, 'y');
        %     ax.YTick = y_tick_values;
        
        str_tit = sprintf('%s - Histogram of Channel #%d', str_tit_1, ch_num);
        title(str_tit)
        
        
        str_save = sprintf('/Results/hist_%s_ch%d.png', str_save_1, ch_num);
        saveas(hh,[pwd str_save]);
        
    end
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