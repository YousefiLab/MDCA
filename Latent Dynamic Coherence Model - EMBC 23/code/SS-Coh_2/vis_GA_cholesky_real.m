function [] = vis_GA_cholesky_real(X, Yk, L, W_iter, MU_iter, w_type, dim, x_tick, x_lim, y_tick, y_lim)

    Y_mapped_chol_init = fun_calc_cholesky(X, Yk, W_iter, L, MU_iter, dim);
    vis_cholesky(Y_mapped_chol_init, w_type,x_tick, x_lim, y_tick, y_lim)
   
end

%% nested functions 
% calc cholesky decomposition 
function Y_mapped_chol = fun_calc_cholesky(x, Yk, W, L, mu, dim)

mm = size(W);
M = length(x);

my = size(Yk);
Y_mapped_chol = zeros(my(2), my(1));



for i=1 : my(1)
    x_k = x(:, i);
    

    D_eval = gen_sigma(x_k, W, dim); % nested fun
    D_eval = inv(D_eval);
    Cov_mat = L*D_eval*L';

    Chol_mat = L*sqrt(D_eval);
    Chol_mat2 = sqrt(D_eval)*L';
%     R = chol(Cov_mat);
    y_k = Yk(i, :);
    yk_remMu = (y_k-mu.');
%     val_1 = yk_remMu*Chol_mat;
%     val_2 = yk_remMu*Chol_mat2;
%     val_3 = Chol_mat*yk_remMu.';
    val_4 = Chol_mat2*yk_remMu.';
    
    Y_mapped_chol(:, i) = val_4;  
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
            lambda_m = exp(x_k(1,1) * w_m(1,1)  +   incongruent_vec(i, 1) *  x_k(2, 1) * w_m(1, 2) + w_bias);
            Sigma_W_x(j, j) = lambda_m;
        end
end

end

% visualize cholesky 
function [] = vis_cholesky(Y_mapped_chol, w_type, x_tick, x_lim, y_tick, y_lim)
nbins = 75;
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
    
    mmm = my(1);
    if mmm> 8
        mmm = 8;
    end
for ch_num_1 =1: mmm
    for ch_num_2 = ch_num_1+1: mmm
        ch_1 = Y_mapped_chol(ch_num_1, :);
        ch_2 = Y_mapped_chol(ch_num_2, :);
        
        min_ = floor(min([ch_1, ch_2]));
        max_ = ceil(max([ch_1, ch_2]));
        
        max_t = max([abs(min_), max_]);
        x_lim_ = [-max_t, max_t];
        x_lim_ = [-8, 8];
        
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