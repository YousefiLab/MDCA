function [mu, L, W, error_tot] = gc_update_mu_L_W(ParamaUpdate, Param, Bayes, Yk, thr_val, count_break, iter, W_sim, update_coef_W, dim)


if nargin < 3
    error('the first 3 inputs "ParamaUpdate, Param, Bayes" should be passed')
end
if nargin == 3 
    thr_val = 0.01;
    count_break = 5;
end
if nargin == 4
    thr_val = 0.01;
    count_break = 5;
end

% extracting data from Bayes
MU_smoother = Bayes.MU_smoother;
SIGMA_smoother = Bayes.SIGMA_smoother;

% extracting from structures
mu_update = ParamaUpdate.mu_update;
W_update = ParamaUpdate.W_update;
L_update = ParamaUpdate.L_update;

dim = Param.dim;
mu = Param.mu;
L = Param.L;
W = Param.W; 
incongruent_vec = Param.incongruent_vec;



error_mu = 0;
error_L = 0;
error_W = 0;

count = 0;
thr = 100;
while thr > thr_val && count < count_break
    count = count+1;
    
    if mu_update == 1
        % estimate mu
        [mu, mu_diff] = gc_estimate_mu(MU_smoother, Yk, mu, W, L, incongruent_vec, dim);
        error_mu(count) = abs(sum(mu_diff(:)));
    end
    
    
    if L_update == 1
        % estimate L
        [L, L_diff] = gc_estimate_L(MU_smoother, SIGMA_smoother, Yk, mu, W, L, incongruent_vec, dim);
        error_L(count) = abs(sum(sum(L_diff)));
    end
    
    if W_update==1
         
        % updating all elements of rows at the same time 
        [W, W_diff]  = gc_estimate_W_2d_1_1(MU_smoother, SIGMA_smoother, Yk, mu, W, L, incongruent_vec,iter, W_sim, update_coef_W, dim);
        
        % updte each element of a row seperately
% %         [W, W_diff]  = gc_estimate_W_2d_1(MU_smoother, SIGMA_smoother, Yk, mu, W, L, incongruent_vec,iter, W_sim, update_coef_W, dim);

        error_W(count) = abs(sum(sum(W_diff)));
    end
    
    error_tot(count) = error_mu(end) + error_L(end) + error_W(end);
    thr = error_tot(count);
    
end

end