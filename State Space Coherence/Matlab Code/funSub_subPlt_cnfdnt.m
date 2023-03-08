function [] = funSub_subPlt_cnfdnt(x_ind, cnfdnc_intvl_filt, MU_filter, cnfdnc_intvl_smoother,MU_smoother,  xk_sample, dim_num, iter_plot)

patch_color = [197, 219, 210; 193, 198, 219]./255;
mu_color = [5, 156, 149; 5, 8, 176]./255; 
x_color = [148, 28, 1; 18, 17, 18]./255;



switch dim_num
    case 1
        x_color = x_color(1, :);
    case 2 
        x_color = x_color(2, :);
end


figure('units','normalized','outerposition',[.25 0.05  .35 .95])
subplot(211)
hold all
patch([x_ind fliplr(x_ind)], [cnfdnc_intvl_filt(1,:) fliplr(cnfdnc_intvl_filt(2,:))], patch_color(1, :), 'LineStyle','none')
plot(MU_filter(1,:), 'Color', mu_color(1, :), 'LineWidth', 2.5);
plot(xk_sample(1, :), 'Color', x_color, 'LineWidth', 3);

xlim([1 x_ind(end)])
str_tit = sprintf('%s -- %d^{st} dimension -- %d Iter ',  'Filter', dim_num , iter_plot);
title(str_tit)
xlabel('trial')

lgd = legend('Confidence Interval', 'Estimated Value', 'X_k simulation');
lgd.FontSize = 5;
hold off


subplot(212)
hold all
patch([x_ind fliplr(x_ind)], [cnfdnc_intvl_smoother(1,:) fliplr(cnfdnc_intvl_smoother(2,:))], patch_color(2, :), 'LineStyle','none')
plot(MU_smoother(1,:), 'Color', mu_color(2, :), 'LineWidth', 2.5);
plot(xk_sample(1, :), 'Color', x_color, 'LineWidth', 3);

xlim([1 x_ind(end)])
str_tit = sprintf('%s -- %d^{st} dimension -- %d Iter ',  'Smoother', dim_num, iter_plot);
title(str_tit)
xlabel('trial')

lgd = legend('Confidence Interval', 'Estimated Value', 'X_k simulation');
lgd.FontSize = 5;

hold off

end