function [W_estimate, dif_W] = gc_estimate_W_2d_2(MU_smoother, SIGMA_smoother, Yk, mu, W, L, incongruent_vec, iter, W_sim)

W_init = W;

m = size(W_init);
M = m(1);

W_estimate = zeros(m);

for m=1: 1
    W_sim_ch = W_sim(m, :);
    
    disp(m)
    
    Wm_init = W_init(m, :);
    Wm1 = Wm_init(1, 1);
    Wm2 = Wm_init(1, 2);
    Wm3 = Wm_init(1, 3);
    
    Um = L(:, m);
    
    x0 = Wm2; % initial guess for solution
    
    x = fminsearch(@(x)f(x, Wm1, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um),x0);
    
    W_estimate(m, :) = [Wm1 x Wm3];
    
    
    if m==1
        func_plot(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, W_sim_ch(2))
    end
end

dif_W = W_estimate - W_init;
end

function obj = f(x, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um)


K = length(MU_smoother);
M = length(mu);


Wm_hat = [x Wm2];
Wm_hat = Wm_hat.';

Wm_biased = Wm3;

Zm = zeros(M);
obj = 0;

for k=1: K
    % p_k - Eq 2.a (UPDATE_L_W.docx)
    y_k = Yk(k, :).';
    p_k = y_k- mu;
    
    mu_s = MU_smoother(:, k);
    sigma_s = squeeze(SIGMA_smoother(k, :, :));
    
    % estimate E[lambda_inv] - Eq 4.c (UPDATE_L_W.docx)
    Exp_term = exp(-((Wm_hat.' * (mu_s-.5*sigma_s*Wm_hat))) - Wm_biased);
    
    obj = obj + Exp_term*p_k'*Um*Um'*p_k;
end

mu_hat = sum(MU_smoother, 2);

obj = obj + Wm_hat.'*mu_hat  +  K*Wm_biased;
% obj = -obj - Wm_hat.'*mu_hat  -  K*Wm_biased;
end


function [] = func_plot(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, w_true)


K = length(MU_smoother);
M = length(mu);


Wm1_vec = Wm1-1: .05 : Wm1 + 1;
Wm1_vec = [Wm1_vec x w_true];
for ww=1: length(Wm1_vec)
    
    Wm1_temp = Wm1_vec(ww);
    
    Wm_hat = [Wm1_temp Wm2];
    Wm_hat = Wm_hat.';
    
    Wm_biased = Wm3;
    
    Zm = zeros(M);
    cost_fun_val = 0;
    
    for k=1: K
        % p_k - Eq 2.a (UPDATE_L_W.docx)
        y_k = Yk(k, :).';
        p_k = y_k- mu;
        
        mu_s = MU_smoother(:, k);
        sigma_s = squeeze(SIGMA_smoother(k, :, :));
        
        % estimate E[lambda_inv] - Eq 4.c (UPDATE_L_W.docx)
        Exp_term = exp(-((Wm_hat.' * (mu_s-.5*sigma_s*Wm_hat))) - Wm_biased);
        
        cost_fun_val = cost_fun_val + Exp_term*p_k'*Um*Um'*p_k;
    end
    
    mu_hat = sum(MU_smoother, 2);
    
    cost_fun_vec(ww) = real(cost_fun_val + Wm_hat.'*mu_hat  +  K*Wm_biased);
    
end

x_val = real(cost_fun_vec(end-1));
w_true_val = real(cost_fun_vec(end));

cost_fun_vec = cost_fun_vec(1:end-2);
Wm1_vec = Wm1_vec(1:end-2);

min_val = min(cost_fun_vec);
ind_min = find(min_val == cost_fun_vec);
Wm1_min = Wm1_vec(ind_min);


str_tit = sprintf('Wm1 -- ROW %d', m);

h=figure;
str_xlabel = sprintf('channel %d', m);
plot(Wm1_vec, cost_fun_vec), xlabel(str_xlabel), ylabel('value of cost function');
title(str_tit)
hold on 
plot(x, x_val,'o', 'MarkerFaceColor', 'r', 'LineWidth', 7)
hold on 
plot(w_true, w_true_val,'*', 'MarkerFaceColor', 'g', 'LineWidth', 10)
legend('cost function', str_xlabel, 'Simulated value')

% filename = sprintf('Wm1_row_%d_iter_%d.jpg', m, iter);
% saveas(h,filename)
% close 
end


