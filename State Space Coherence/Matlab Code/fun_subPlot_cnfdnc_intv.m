function [] = fun_subPlot_cnfdnc_intv(xk_sample, MU_filter, SIGMA_filter, MU_smoother, SIGMA_smoother, dim, iter_plot)

cnfdnc_intvl_filt = get_cnfdnt_intvl(MU_filter, SIGMA_filter, dim);
cnfdnc_intvl_smoother = get_cnfdnt_intvl(MU_smoother, SIGMA_smoother, dim);


x_ind = 1: length(cnfdnc_intvl_filt);

switch dim
    case 1
        cnfdnc_intvl_filt_d1 = cnfdnc_intvl_filt(1:2, :);
        cnfdnc_intvl_smoother_d1 = cnfdnc_intvl_smoother(1:2, :);
        
        % call plot function
        funSub_subPlt_cnfdnt(x_ind, cnfdnc_intvl_filt_d1, MU_filter, cnfdnc_intvl_smoother_d1, MU_smoother,  xk_sample, 1, iter_plot)
        
    case 2
        cnfdnc_intvl_filt_d1 = cnfdnc_intvl_filt(1:2, :);
        cnfdnc_intvl_filt_d2 = cnfdnc_intvl_filt(3:4, :);
        cnfdnc_intvl_smoother_d1 = cnfdnc_intvl_smoother(1:2, :);
        cnfdnc_intvl_smoother_d2 = cnfdnc_intvl_smoother(3:4, :);
        
        % call plot function
        funSub_subPlt_cnfdnt(x_ind(), cnfdnc_intvl_filt_d1, MU_filter(1, :), cnfdnc_intvl_smoother_d1, MU_smoother(1, :),  xk_sample(1, :), 1, iter_plot)
        funSub_subPlt_cnfdnt(x_ind, cnfdnc_intvl_filt_d2, MU_filter(2, :), cnfdnc_intvl_smoother_d2, MU_smoother(2, :),  xk_sample(2, :), 2, iter_plot)
        
       
end


end

