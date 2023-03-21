function [xk_sample, Yk_sample] = SIM_samp_x_yk_switch(cnfg_sim, switch_num)
dim = cnfg_sim.dim;
x0 = cnfg_sim.x0;
F = cnfg_sim.F;
D = cnfg_sim.D;
Q = cnfg_sim.Q;
incongruent_vec = cnfg_sim.incongruent_vec;
L_ = cnfg_sim.L;
L = L_{switch_num};
W = cnfg_sim.W;
mu = cnfg_sim.mu;

switch  dim
    case 1
        x0 = x0(1,1);
        F = F(1,1);
        D = D(1,1);
        Q = Q(1, 1);
    case 2
        
end

M = length(incongruent_vec);

%% generate samples of x
xk_sample = zeros(dim, M+1);
MU_x = zeros(dim, M);

xk_sample(:, 1) = x0;
for i=1: M
    
    MU_x(:, i) = F*xk_sample(:, i) + D; % mu of normal distribution
    
    % calculate normal pdf for elemants of x with a given mu (x_0) and var (s_v)
    switch dim
        %1d
        case 1
            xk_sample(:, i+1) = normrnd(MU_x(:, i) ,Q);
            
            %2d
        case 2
            xk_sample(:, i+1) = mvnrnd(MU_x(:, i), Q);
    end
end
xk_sample(:, 1) = [];

% step = 69;
% xk_sample_2(1, 1:280)=0;
% xk_sample_2(1, 281:350)= 0:1/step:1;
% xk_sample_2(1, 300:700)=1;
% xk_sample_2(1, 701:770)= 1:-1/step:0;
% xk_sample_2(1, 721:1000)=0;
% 
% xk_sample = xk_sample+ xk_sample_2;

%% generate sample of SIGMA
mm = size(W);
mu_sigma = zeros(mm(1), 1);
Sigma_sample = zeros(mm(1), M);
Yk_sample = zeros(M, mm(1));
for i=1 : M
    x_temp = xk_sample(:, i);
    
    switch dim
        
        %1 d
        case 1
            Sigma_W_x = zeros(mm(1));
            for j=1: mm(1)
                w_m = W(j, 1);
                w_bias = W(j,3);
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
    
    
    Sigma_sample(:, i) = mvnrnd(mu_sigma, Sigma_W_x, 1);
    Yk_sample(i, :) = mu + L*Sigma_sample(:, i);
    
    
end

end