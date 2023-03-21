function [Param_Update, error_tot] = gc_parameter_update(Param, ParamaUpdate, Bayes, Yk, iter, W_sim, update_coef_W)
% GC_PARAMETER_UPDATE this function is used to estimate parameters


if nargin < 3
    error('the first 3 inputs "ParamaUpdate, Param, Bayes" should be passed')
end


% strcutur of Param
F = Param.F;
Q = Param.Q;
dim = Param.dim;

% update of don't update info
sigma_update = ParamaUpdate.sigma_update;
mu_update = ParamaUpdate.mu_update;
W_update = ParamaUpdate.W_update;
L_update = ParamaUpdate.L_update;

Param_Update = Param;

if sigma_update == 1
    % estimate sv2
    [Q_estimate, error_tot] = gc_update_sigma(Bayes, F, Q, dim);
    Param_Update.Q = Q_estimate;
end

% update mu, L and am & bm
if (mu_update+W_update+L_update) >= 1
    thr_val = 1e-6;
    count_break = 1;
    [mu_estimate, L_estimate, W_estimate, error_tot] = gc_update_mu_L_W(ParamaUpdate, Param, Bayes, Yk, thr_val, count_break, iter, W_sim, update_coef_W, dim);
    Param_Update.mu = mu_estimate;
    Param_Update.W = W_estimate;
    Param_Update.L = L_estimate;
end

end