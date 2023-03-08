function Yk_sample = gen_yk_of_xk(xk, W, L, mu)

mm = size(W);
M = length(xk);
mu_sigma = zeros(mm(1), 1);
Sigma_sample = zeros(mm(1), M);
Yk_sample = zeros(M, mm(1));
dim = mm(2)-1;

for i=1 : M
    x_temp = xk(:, i);
    
    switch dim
        
        %1 d
        case 1
            Sigma_W_x = zeros(mm(1));
            for j=1: mm(1)
                w_m = W(j, 1);
                w_bias = W(j,2);
                lambda_m = exp(x_temp(1,1) * w_m(1,1)+ w_bias);
                Sigma_W_x(j, j) = lambda_m;
            end
            
            
            % 2d
        case 2
            Sigma_W_x = zeros(mm(1));
            for j=1: mm(1)
                w_m = W(j, 1:2);
                w_bias = W(j,3);
                lambda_m = exp(x_temp(1,1) * w_m(1,1)  +   incongruent_vec(i, 1) *  x_temp(2, 1) * w_m(1, 2) + w_bias);
                Sigma_W_x(j, j) = lambda_m;
            end
            
            
    end
    
    
    Sigma_sample(:, i) = mvnrnd(mu_sigma, Sigma_W_x);
    Yk_sample(i, :) = mu + L*Sigma_sample(:, i);
    
    
end
end