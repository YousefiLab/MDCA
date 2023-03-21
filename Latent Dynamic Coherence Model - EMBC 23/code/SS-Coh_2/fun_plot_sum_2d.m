function []= fun_plot_sum_2d(xk_sample,MU_filter, MU_smoother, str_tit)

figure, 
hold all 
plot(1: length(xk_sample), sum(xk_sample), 'LineWidth', 3)
plot(1: length(xk_sample), sum(MU_filter), 'LineWidth', 3)
plot(1: length(xk_sample), sum(MU_smoother), 'LineWidth', 3)
hold off
legend('simulated x_k', 'filter estimatoin', 'smoother estimation')
xlim([1 length(xk_sample)])


title(str_tit)



end