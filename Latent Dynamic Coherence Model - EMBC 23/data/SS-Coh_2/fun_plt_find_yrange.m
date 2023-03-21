function str_min_max = fun_plt_find_yrange(bayes_1, bayes_2, xk_sample, dim, using_sim_data)

if using_sim_data==1
    cnfdnc_intvl_1 = cnfdnt_fun(bayes_1, dim);
    cnfdnc_intvl_2 = cnfdnt_fun(bayes_2, dim);
    
    str_min_max = [];
    switch dim
        case 1
            % for first bayes
            cnfdnc_intvl_1_d1 = cnfdnc_intvl_1(1:2, :);
            cnfdnc_intvl_2_d1 = cnfdnc_intvl_2(1:2, :);
            
            smoother_and_x_d1 = [cnfdnc_intvl_1_d1; cnfdnc_intvl_2_d1; xk_sample];
            
            max_1_d1 = max(smoother_and_x_d1(:));
            min_1_d1 = min(smoother_and_x_d1(:));
            
            str_min_max.min_1_d1 = min_1_d1;
            str_min_max.max_1_d1 = max_1_d1;
        case 2
            % for first bayes
            cnfdnc_intvl_1_d1 = cnfdnc_intvl_1(1:2, :);
            cnfdnc_intvl_1_d2 = cnfdnc_intvl_1(3:4, :);
            
            % for second bayes
            cnfdnc_intvl_2_d1 = cnfdnc_intvl_2(1:2, :);
            cnfdnc_intvl_2_d2 = cnfdnc_intvl_2(3:4, :);
            
            smoother_and_x_d1 = [cnfdnc_intvl_1_d1; cnfdnc_intvl_2_d1; xk_sample(1,:)];
            smoother_and_x_d2 = [cnfdnc_intvl_1_d1; cnfdnc_intvl_2_d1; xk_sample(2,:)];
            
            max_1_d1 = max(smoother_and_x_d1(:));
            min_1_d1 = min(smoother_and_x_d1(:));
            
            max_1_d2 = max(smoother_and_x_d2(:));
            min_1_d2 = min(smoother_and_x_d2(:));
            
            str_min_max.min_1_d1 = min_1_d1;
            str_min_max.max_1_d1 = max_1_d1;
            str_min_max.min_1_d2 = min_1_d2;
            str_min_max.max_1_d2 = max_1_d2;
    end
    
elseif using_sim_data==0
    cnfdnc_intvl_1 = cnfdnt_fun(bayes_1, dim);
    cnfdnc_intvl_2 = cnfdnt_fun(bayes_2, dim);
    
    str_min_max = [];
    switch dim
        case 1
            % for first bayes
            cnfdnc_intvl_1_d1 = cnfdnc_intvl_1(1:2, :);
            cnfdnc_intvl_2_d1 = cnfdnc_intvl_2(1:2, :);
            
            smoother_and_x_d1 = [cnfdnc_intvl_1_d1; cnfdnc_intvl_2_d1];
            
            max_1_d1 = max(smoother_and_x_d1(:));
            min_1_d1 = min(smoother_and_x_d1(:));
            
            str_min_max.min_1_d1 = min_1_d1;
            str_min_max.max_1_d1 = max_1_d1;
        case 2
            % for first bayes
            cnfdnc_intvl_1_d1 = cnfdnc_intvl_1(1:2, :);
            cnfdnc_intvl_1_d2 = cnfdnc_intvl_1(3:4, :);
            
            % for second bayes
            cnfdnc_intvl_2_d1 = cnfdnc_intvl_2(1:2, :);
            cnfdnc_intvl_2_d2 = cnfdnc_intvl_2(3:4, :);
            
            smoother_and_x_d1 = [cnfdnc_intvl_1_d1; cnfdnc_intvl_2_d1];
            smoother_and_x_d2 = [cnfdnc_intvl_1_d1; cnfdnc_intvl_2_d1];
            
            max_1_d1 = max(smoother_and_x_d1(:));
            min_1_d1 = min(smoother_and_x_d1(:));
            
            max_1_d2 = max(smoother_and_x_d2(:));
            min_1_d2 = min(smoother_and_x_d2(:));
            
            str_min_max.min_1_d1 = min_1_d1;
            str_min_max.max_1_d1 = max_1_d1;
            str_min_max.min_1_d2 = min_1_d2;
            str_min_max.max_1_d2 = max_1_d2;
    end
    
end


end

function cnfdnc_intvl = cnfdnt_fun(bayes_strct, dim)

mu_estimation = bayes_strct.MU_smoother;
sigma_estimation = bayes_strct.SIGMA_smoother;

m = size(mu_estimation);

cnfdnc_intvl = zeros(2*m(1), m(2));
% 1st row: up confidence interval of 1d
% 2nd row: down confidence interval of 1d
% 3rd row: up confidence interval of 2d
% 4th row: down confidence interval of 2d

for i=1: m(1) % number of dimensiona (2)
    
    for j=1: m(2) % number of points (127)
        mu_temp_dim = mu_estimation(i, j);
        sigma_temp_dim = squeeze(sigma_estimation(j, :, :));
        std_temp = sqrt(sigma_temp_dim(i, i));
        
        up_lim_temp = mu_temp_dim + 2*std_temp;
        down_lim_temp = mu_temp_dim - 2*std_temp;
        
        cnfdnc_intvl( (i-1)*dim + 1, j ) = up_lim_temp;
        cnfdnc_intvl( (i-1)*dim + 2, j ) = down_lim_temp;
    end
end




end
