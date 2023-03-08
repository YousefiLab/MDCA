function [mu_update, diff_mu] = gc_estimate_mu(MU_smoother, Y_k, mu, W, L, incongruent_vec, dim)


mu_init = mu;

% get expected of D^-1
E_Sigma_inv = gc_expectation_D_inv(MU_smoother, W, incongruent_vec, dim);


m = size(Y_k);
% m(1) = number of steps (K)
% m(2) = number of channels

sum_num = zeros(m(2), 1);
sum_denum = zeros(m(2), m(2));

for k=1 : m(1)
    E_inv_k = squeeze(E_Sigma_inv(k, :,:));
    Yk_k = Y_k(k,:);
    
    sum_num = sum_num + (L*E_inv_k*L')*Yk_k.';
    sum_denum = sum_denum + L*E_inv_k*L';
    
end

mu_update = inv(sum_denum)*sum_num;

diff_mu = mu_update - mu_init;

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate expectation of D^-1
function [E_Sigma_W_x_inv] = gc_expectation_D_inv(MU_smoother, W, incongruent_vec, dim)

K = length(MU_smoother);
sigma_dim = length(W);
Sigma_inv = zeros(K, sigma_dim, sigma_dim);

for k=1 : K
    x_k = MU_smoother(:, k);
    
    
    switch dim
        
        %1 d
        case 1
            Sigma_W_x_inv = zeros(sigma_dim);
            
            for j=1: sigma_dim
                w_m = W(j, 1);
                w_bias = W(1,2);
                lambda_m = exp(x_k(1,1) * w_m(1,1)+ w_bias);
                Sigma_W_x_inv(j, j) = inv(lambda_m);
            end
            
            
            % 2d
        case 2
            Sigma_W_x_inv = zeros(sigma_dim);
            
            for j=1: sigma_dim
                w_m = W(j, 1:2);
                w_bias = W(1,3);
                lambda_m = exp(x_k(1,1) * w_m(1,1)  +   incongruent_vec(k, 1) *  x_k(2, 1) * w_m(1, 2) + w_bias);
                Sigma_W_x_inv(j, j) = inv(lambda_m);
            end
            
    end
    
    Sigma_inv(k, :, :) = Sigma_W_x_inv;
    
end
E_Sigma_W_x_inv = Sigma_inv;
end