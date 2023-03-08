function [MU_filter, SIGMA_filter, MU_oneStep, SIGMA_oneStep] = gc_bayes_filter_oneStep(Y_k, L, mu, incongruent_vec, W, mu_1_0, sigma_1_0, F, D, Q, dim)
% GC_BAYES_FILTER this function estimate the posterior probability of the
% state by estimating the Likelihood and One-Step Prediction distributions
% (for further information please take a look at part 1-5 of supplemntary)
%
% INPUTS
% f_xx      : state transient matrix
% p_xx_pre  : distribution of previous posterior
% y_k       : the current observation
% x         : an interval of state variables
% L         : given eigenvector
% mu        : given mu
% incongruent_vec : 1 trail is incngruent, 0 trial is congruent
% W               : matrix of weight
% mu_1_0          : initial value for mu of filtre used in GA
% sigma_1_0       : initial value for sigma of filter used in GA


switch dim
    case 1
        mu_1_0 = mu_1_0(1,1);
        sigma_1_0 = sigma_1_0(1,1);
        F = F(1,1);
        D = D(1,1);
        Q = Q(1,1);
    case 2
        
end


m = size(Y_k);
% m(1) = K (number of states)
% m(2) = number of channels
K = m(1);






MU_filter = zeros(dim, K);
SIGMA_filter = zeros(K, dim, dim);
MU_oneStep = zeros(dim, K);
SIGMA_oneStep = zeros(K, dim, dim);

MU_oneStep(:, 1) = mu_1_0;
SIGMA_oneStep(1, :, :) = sigma_1_0;


% % use squeeze(SIGMA_FILTER(1, :, :)) to convert 3d into 2d

% estimation of posterior probablity of x based BAYES FILTER
for i=1:K
%     disp(i)
    y_k = Y_k(i,:).';  % current observation
    

    mu_k_1 = MU_oneStep(:, i);
    sigma_k_1 = squeeze(SIGMA_oneStep(i, :, :));

    % calculate U
%     Ly = L*(y_k-mu);
    Ly = conj(L')*(y_k-mu);
    U = Ly.*conj(Ly);
   
    
    m = size(Ly);
    

    
    %% estimation of sigma 
    switch dim
        case 1
           sigma_W_lambda_U_m = zeros(dim);
            for j=1 : m(1)
                w_m = W(j, 1);
                w_bias = W(j,2);
                lambda_m = exp(mu_k_1(1,1) * w_m(1,1) + w_bias);
                sigma_W_lambda_U_m = sigma_W_lambda_U_m + w_m.^2 * inv(lambda_m) * U(j, 1);
            end 
            
            
        case 2
            sigma_W_lambda_U_m = zeros(dim);
            for j=1 : m(1)
                w_m = W(j, 1:2);
                w_bias = W(j,3);
                lambda_m = exp(mu_k_1(1,1) * w_m(1,1) + incongruent_vec(i, 1) *  mu_k_1(2, 1) * w_m(1, 2) + w_bias);
                sigma_W_lambda_U_m = sigma_W_lambda_U_m + w_m.' * w_m * inv(lambda_m) * U(j, 1);
            end
    end
      
    
    SIGMA_FILTER_inv = inv(sigma_k_1) + sigma_W_lambda_U_m;
    sigma_filter = inv(SIGMA_FILTER_inv);
    SIGMA_filter(i, :, :) = sigma_filter;
    
    
    
    %% estimation of mu
    
    
        switch dim
            %1 d
        case 1
            sigma_task = zeros(m(1), m(1));
            for j=1: m(1)
                w_m = W(j, 1);
                w_bias = W(j,2);
                lambda_m = exp(mu_k_1(1,1)*w_m(1,1) + w_bias);
                sigma_task(j, j) = lambda_m;
            end
            W_temp = W(:, 1);
            
            
            % 2d
        case 2
            sigma_task = zeros(m(1), m(1));
            for j=1: m(1)
                w_m = W(j, 1:2);
                w_bias = W(j,3);
                lambda_m = exp(mu_k_1(1,1)*w_m(1,1) + incongruent_vec(i, 1)* mu_k_1(2, 1)*w_m(1, 2) + w_bias);
                sigma_task(j, j) = lambda_m;
            end
            W_temp = W(:, 1:2);
            
        end
    
    
    mu_filter = mu_k_1 + sigma_filter * (W_temp.' * inv(sigma_task) * U - W_temp.' * ones(m(1),1));
    MU_filter(:, i) =  mu_filter;
    
    
    mu_k_1 = MU_filter(:, i);
    sigma_k_1 = squeeze(SIGMA_filter(i, :, :));
    
    %% update oneStep 
    MU_oneStep(:, i+1) = F*mu_k_1 + D;
    SIGMA_oneStep(i+1, :, :) = F*sigma_k_1*F.' + Q;
end
SIGMA_oneStep(i+1, :, :) = [];


end

