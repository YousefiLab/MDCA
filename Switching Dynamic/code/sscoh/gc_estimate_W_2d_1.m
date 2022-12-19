function [W_estimate, dif_W] = gc_estimate_W_2d_1(MU_smoother, SIGMA_smoother, Yk, mu, W, L, incongruent_vec, iter, W_sim, update_coef_W, dim)

W_init = W;

m = size(W_init);
ch_num = m(1);
Wm_num = m(2);

W_estimate = W;

for m=1: ch_num
    W_sim_ch = W_sim(m, :);
    
    Um = L(:, m);
    
    for ww=1: Wm_num
        %     disp(m)
        
        % check we want to update or not
        updata_val = update_coef_W(m, ww);
        
        if updata_val==1
            
            %% for 1d 
            if dim==1
                Wm_init = W_init(m, :);
                Wm1 = Wm_init(1, 1);
                Wm3 = Wm_init(1, 2);
                
                switch ww
                    
                    case 1 %x1
                        
                        x0 = Wm1; % initial guess for solution
                        
                        % call find min function
                        x = fminsearch(@(x)f_update_wm1_1d(x, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um),x0);
                        
                        % update W_estimate
                        W_estimate(m, ww) = x;
                        
                        % plot cost function
                        % % %                     func_plot_wm1(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, W_sim_ch(1))
                        
                        
                    case 2 %bias
                        
                        x0 = Wm3; % initial guess for solution
                        
                        % call find min function
                        x = fminsearch(@(x)f_update_wm3_1d(x, Wm1, MU_smoother, SIGMA_smoother, Yk, mu, Um),x0);
                        
                        % update W_estimate
                        W_estimate(m, ww) = x;
                        
                        % plot cost function
                        % % %                     func_plot_wm3(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, W_sim_ch(3))
                        
                        
                end
            end
            
            %% for 2d 
            if dim==2
                
                Wm_init = W_init(m, :);
                Wm1 = Wm_init(1, 1);
                Wm2 = Wm_init(1, 2);
                Wm3 = Wm_init(1, 3);
                
                
                switch ww
                    case 1 %  x1
                        
                        x0 = Wm1; % initial guess for solution
                        
                        % call find min function
                        x = fminsearch(@(x)f_update_wm1_2d(x, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um),x0);
                        
                        % update W_estimate
                        W_estimate(m, ww) = x;
                        
                        % plot cost function
                        % % %                     func_plot_wm1(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, W_sim_ch(1))
                        
                        
                    case 2 % x2 
                        
                        x0 = Wm2; % initial guess for solution
                        
                        % call find min function
                        x = fminsearch(@(x)f_update_wm2_2d(x, Wm1, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um),x0);
                        
                        % update W_estimate
                        W_estimate(m, ww) = x;
                        
                        %                     if m==1
                        % plot cost function
                        % %                         func_plot_wm2(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, W_sim_ch(2))
                        %                     end
                        
                    case 3 % bias
                        
                        x0 = Wm3; % initial guess for solution
                        
                        % call find min function
                        x = fminsearch(@(x)f_update_wm3_2d(x, Wm1, Wm2, MU_smoother, SIGMA_smoother, Yk, mu, Um),x0);
                        
                        % update W_estimate
                        W_estimate(m, ww) = x;
                        
                        % plot cost function
                        % % %                     func_plot_wm3(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, W_sim_ch(3))
                end
            end
            
            
        end
        
    end
    
    
end

dif_W = W_estimate - W_init;
end

% nested function for w of x1 1d
function obj = f_update_wm1_1d(x, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um)


K = length(MU_smoother);
M = length(mu);


Wm_hat = x;
Wm_hat = Wm_hat.';

Wm_biased = Wm3;

Zm = zeros(M);
obj = 0;

for k=1: K
    % p_k - Eq 2.a (UPDATE_L_W.docx)
    y_k = Yk(k, :).';
    p_k = y_k- mu;
    
    mu_s = MU_smoother(k);
    sigma_s = SIGMA_smoother(k);
    
    % estimate E[lambda_inv] - Eq 4.c (UPDATE_L_W.docx)
    Exp_term = exp(-((Wm_hat.' * (mu_s-.5*sigma_s*Wm_hat))) - Wm_biased);
    
    obj = obj + Exp_term*p_k'*Um*Um'*p_k;
end

mu_hat = sum(MU_smoother, 2);

obj = obj + Wm_hat.'*mu_hat  +  K*Wm_biased;
% obj = -obj - Wm_hat.'*mu_hat  -  K*Wm_biased;
end

% nested function for w bias 1d
function obj = f_update_wm3_1d(x, Wm1, MU_smoother, SIGMA_smoother, Yk, mu, Um)


K = length(MU_smoother);
M = length(mu);


Wm_hat = [Wm1];
Wm_hat = Wm_hat.';


obj = 0;

for k=1: K
    % p_k - Eq 2.a (UPDATE_L_W.docx)
    y_k = Yk(k, :).';
    p_k = y_k- mu;
    
    mu_s = MU_smoother(k);
    sigma_s = SIGMA_smoother(k);
    
    % estimate E[lambda_inv] - Eq 4.c (UPDATE_L_W.docx)
    Exp_term = exp(-((Wm_hat.' * (mu_s-.5*sigma_s*Wm_hat))) - x);
    
    obj = obj + Exp_term*p_k'*Um*Um'*p_k;
end

mu_hat = sum(MU_smoother, 2);

obj = obj + Wm_hat.'*mu_hat  +  K*x;
% obj = -obj;
end

% nested function for w of x1 2d
function obj = f_update_wm1_2d(x, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um)


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

% nested function for w of x2 2d
function obj = f_update_wm2_2d(x, Wm1, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um)


K = length(MU_smoother);
M = length(mu);


Wm_hat = [Wm1 x];
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
% obj = -obj;
end

% nested function for w bias 2d
function obj = f_update_wm3_2d(x, Wm1, Wm2, MU_smoother, SIGMA_smoother, Yk, mu, Um)


K = length(MU_smoother);
M = length(mu);


Wm_hat = [Wm1 Wm2];
Wm_hat = Wm_hat.';


obj = 0;

for k=1: K
    % p_k - Eq 2.a (UPDATE_L_W.docx)
    y_k = Yk(k, :).';
    p_k = y_k- mu;
    
    mu_s = MU_smoother(:, k);
    sigma_s = squeeze(SIGMA_smoother(k, :, :));
    
    % estimate E[lambda_inv] - Eq 4.c (UPDATE_L_W.docx)
    Exp_term = exp(-((Wm_hat.' * (mu_s-.5*sigma_s*Wm_hat))) - x);
    
    obj = obj + Exp_term*p_k'*Um*Um'*p_k;
end

mu_hat = sum(MU_smoother, 2);

obj = obj + Wm_hat.'*mu_hat  +  K*x;
% obj = -obj;
end

function [] = func_plot_wm1(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, w_true)


K = length(MU_smoother);
M = length(mu);

itv_plt = 2;
Wm1_vec = Wm1-itv_plt: .1 : Wm1 + itv_plt;
Wm1_vec = [Wm1_vec x w_true];
for ww=1: length(Wm1_vec)
    
    Wm1_temp = Wm1_vec(ww);
    
    Wm_hat = [Wm1_temp Wm2];
    Wm_hat = Wm_hat.';
    
    Wm_biased = Wm3;
    
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
    
    cost_fun_val = real(cost_fun_val + Wm_hat.'*mu_hat  +  K*Wm_biased);
    cost_fun_vec(ww) = cost_fun_val;
end

x_val = real(cost_fun_vec(end-1));
w_true_val = real(cost_fun_vec(end));

cost_fun_vec = cost_fun_vec(1:end-2);
Wm1_vec = Wm1_vec(1:end-2);


str_tit = sprintf('Wm1 -- ch %d -- iter %d', m, iter);

h=figure;
str_xlabel = sprintf('channel %d', m);
plot(Wm1_vec, cost_fun_vec), xlabel(str_xlabel), ylabel('value of cost function');
title(str_tit)
hold on
plot(x, x_val,'o', 'MarkerFaceColor', 'r', 'LineWidth', 7)
hold on
plot(w_true, w_true_val,'*', 'MarkerFaceColor', 'g', 'LineWidth', 10)
legend('cost function', str_xlabel, 'Simulated value')

% filename = sprintf('Wm1_ch_%d_iter_%d.jpg', m, iter);
% saveas(h,filename)
% close
end


function [] = func_plot_wm2(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, w_true)


K = length(MU_smoother);
M = length(mu);

itv_plt = 2;
Wm2_vec = Wm2-itv_plt: .1 : Wm2 + itv_plt;
Wm2_vec = [Wm2_vec x w_true];
for ww=1: length(Wm2_vec)
    
    Wm2_temp = Wm2_vec(ww);
    
    Wm_hat = [Wm1 Wm2_temp];
    Wm_hat = Wm_hat.';
    
    Wm_biased = Wm3;
    
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
    
    cost_fun_val = real((cost_fun_val + Wm_hat.'*mu_hat  +  K*Wm_biased));
    cost_fun_vec(ww) = cost_fun_val;
end

x_val = real(cost_fun_vec(end-1));
w_true_val = real(cost_fun_vec(end));

cost_fun_vec = cost_fun_vec(1:end-2);
Wm2_vec = Wm2_vec(1:end-2);



str_tit = sprintf('Wm2 -- ch %d -- iter %d', m, iter);

h=figure;
str_xlabel = sprintf('channel %d', m);
plot(Wm2_vec, cost_fun_vec), xlabel(str_xlabel), ylabel('value of cost function');
title(str_tit)
hold on
plot(x, x_val,'o', 'MarkerFaceColor', 'r', 'LineWidth', 7)
hold on
plot(w_true, w_true_val,'*', 'MarkerFaceColor', 'g', 'LineWidth', 10)
legend('cost function', str_xlabel, 'Simulated value')

filename = sprintf('Wm2_ch_%d_iter_%d.jpg', m, iter);
saveas(h,filename)
close
end

function [] = func_plot_wm3(Wm1, Wm2, Wm3, MU_smoother, SIGMA_smoother, Yk, mu, Um, x, m, iter, w_true)


K = length(MU_smoother);
M = length(mu);

itv_plt = 2;
Wm3_vec = Wm3-itv_plt: .1 : Wm3 + itv_plt;
Wm3_vec = [Wm3_vec x w_true];
for ww=1: length(Wm3_vec)
    
    Wm3_temp = Wm3_vec(ww);
    
    Wm_hat = [Wm1 Wm2];
    Wm_hat = Wm_hat.';
    
    Wm_biased = Wm3_temp;
    
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
    
    cost_fun_val = real((cost_fun_val + Wm_hat.'*mu_hat  +  K*Wm_biased));
    cost_fun_vec(ww) = cost_fun_val;
    
end

x_val = real(cost_fun_vec(end-1));
w_true_val = real(cost_fun_vec(end));

cost_fun_vec = cost_fun_vec(1:end-2);
Wm3_vec = Wm3_vec(1:end-2);


str_tit = sprintf('Wm3 -- ch %d -- iter %d', m, iter);

h=figure;
str_xlabel = sprintf('channel %d', m);
plot(Wm3_vec, cost_fun_vec), xlabel(str_xlabel), ylabel('value of cost function');
title(str_tit)
hold on
plot(x, x_val,'o', 'MarkerFaceColor', 'r', 'LineWidth', 7)
hold on
plot(w_true, w_true_val,'*', 'MarkerFaceColor', 'g', 'LineWidth', 10)
legend('cost function', str_xlabel, 'Simulated value')

filename = sprintf('Wm3_ch_%d_iter_%d.jpg', m, iter);
saveas(h,filename)
close
end


