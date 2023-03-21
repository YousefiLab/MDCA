function [cnfdnc_intvl] = get_cnfdnt_intvl(mu_estimation, sigma_estimation, dim)

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