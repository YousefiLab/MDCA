function [sigma_update, sigma_error] = gc_update_sigma_switch(Bayes, F, Q, dim, mat_switch_active_ind)


MU_smoother = Bayes.MU_smoother;
SIGMA_smoother = Bayes.SIGMA_smoother;
G_k = Bayes.G_k;


%% extract training data based on switch 
vec_active = mat_switch_active_ind(1, :);
ind_active = find(vec_active==1);
cntr = 0;
for i_ =1: length(ind_active)
    seg_ = ind_active(i_);
    ind_l = mat_switch_active_ind(2, seg_);
    ind_u = mat_switch_active_ind(3, seg_);
    cntr = cntr +ind_u-ind_l +1;
end

MU_smoother_n = zeros(1, cntr); 
SIGMA_smoother_n = zeros(cntr, 1);
G_k_n = zeros(cntr, 1);


cntr = 0;
for i_ =1: length(ind_active)
    seg_ = ind_active(i_);
    ind_l = mat_switch_active_ind(2, seg_);
    ind_u = mat_switch_active_ind(3, seg_);
    
    l_n = ind_u-ind_l +1;

    MU_smoother_n(1, cntr+1: cntr+l_n) = MU_smoother(1, ind_l: ind_u);
    SIGMA_smoother_n(cntr+1: cntr+l_n, 1) = SIGMA_smoother(ind_l: ind_u, 1);
    G_k_n(cntr+1: cntr+l_n, 1) = G_k(ind_l: ind_u, 1);
    
    cntr = cntr+l_n;
end

%% 

switch dim
    case 1
        F = F(1,1);
        Q = Q(1,1);
    case 2
        
end

% calculate Expectation 
Q_update = zeros(dim);
K = length(SIGMA_smoother_n);

for k=2: K
    
    sigma_smoother_k = squeeze(SIGMA_smoother_n(k, :, :));
    sigma_smoother_k_1 = squeeze(SIGMA_smoother_n(k-1, :, :));
    
    mu_smoother_k = MU_smoother_n(:,k);
    mu_smoother_k_1 = MU_smoother_n(:,k-1);
    
    % from Bayesian filtering and smoothing - sarkka book 
    % this factor has been calculated in smoother step
    Gk_k_1 = squeeze(G_k_n(k-1, :, :)); 
    
    E_xk_xk = gc_expected_x2(sigma_smoother_k , mu_smoother_k);
    E_xk_1_xk_1 = gc_expected_x2( sigma_smoother_k_1, mu_smoother_k_1);
    [expected_joint_xk_xk_1, expected_joint_xk_1_xk] = gc_joint_expectation(mu_smoother_k, mu_smoother_k_1, sigma_smoother_k, Gk_k_1);
    
    Q_update = Q_update +  E_xk_xk + E_xk_1_xk_1 - (expected_joint_xk_xk_1 + expected_joint_xk_1_xk);

end


sigma =  Q_update./(K-1);

% remove off diagnoal terms 
sigma_update = sigma;
% sigma_update = diag(diag(sigma));

sigma_error = abs(sum(sum(sigma_update - Q)));
end








%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calculate joint expectation 
function [expected_joint_k_k_1, expected_joint_k_1_k] = gc_joint_expectation(mu_k, mu_k_1, sigma_smoother_k, Gk_k_1)

sigma_joint_k_k_1 = sigma_smoother_k*Gk_k_1';
sigma_joint_k_1_k = Gk_k_1*sigma_smoother_k;

expected_joint_k_k_1 = mu_k*mu_k_1' + sigma_joint_k_k_1;
expected_joint_k_1_k = mu_k_1*mu_k' + sigma_joint_k_1_k;
end

%% calculate expectation xk.^2
function [expected_x_2] = gc_expected_x2(sigma, mu)

expected_x_2 = sigma + mu*mu.';

end





