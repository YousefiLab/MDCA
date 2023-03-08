function [MU_smoother, SIGMA_smoother, G_k]= gc_bayes_smoother(MU_filter, SIGMA_filter, MU_oneStep, SIGMA_oneStep, F, Q, dim )
% GC_BAYES_SMOOTHER used to estimate the posterior probability based Bayes

switch dim
    case 1
        F = F(1,1);
        
    case 2
        
end


m_mu = size(MU_filter);
m_sigma = size(SIGMA_filter);
K = m_mu(2);

MU_smoother = zeros(m_mu(1), m_mu(2));
MU_smoother(:, K)= MU_filter(:, K);

switch dim
    % d1
    case 1 
        SIGMA_smoother = zeros(m_sigma(1), m_sigma(2));
        G_k = zeros(m_sigma(1)-1, m_sigma(2));
        SIGMA_smoother(K, 1)= SIGMA_filter(K, 1);
        
        % d2
    case 2
        SIGMA_smoother = zeros(m_sigma(1), m_sigma(2), m_sigma(3));
        G_k = zeros(m_sigma(1)-1, m_sigma(2), m_sigma(3));
        
        SIGMA_smoother(K, :, :)= squeeze(SIGMA_filter(K, :, :));
        
end



for i=K-1:-1: 1
    mu_filter_curr = MU_filter(:, i);
    sigma_filter_curr = squeeze(SIGMA_filter(i, :, :));
    
    mu_oneStep_next = MU_oneStep(:, i+1);
    sigma_oneStep_next = squeeze(SIGMA_oneStep(i+1, :, :));
    
    mu_smoother_next = MU_smoother(:, i+1);
    sigma_smoother_next= squeeze(SIGMA_smoother(i+1, :, :));
    
    Gk = sigma_filter_curr*F.'*inv(sigma_oneStep_next);
    G_k(i, :, :) = Gk;
    
    mu_smoother_curr = mu_filter_curr + Gk*(mu_smoother_next - mu_oneStep_next);
    simga_smoother_curr = sigma_filter_curr + Gk* (sigma_smoother_next- sigma_oneStep_next)*Gk';
    
    MU_smoother(:, i)= mu_smoother_curr;
    SIGMA_smoother(i, :, :)= simga_smoother_curr;
    
    
end

end


