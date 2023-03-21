function [likeLihood_val] = gc_likelihood_val(param, Bayes)

incongruent_vec = param.incongruent_vec;
mu = param.mu;
L = param.L;
W = param.W;
F = param.F;
D = param.D;
dim = param.dim;
Q = param.Q;

mu_smoother = Bayes.MU_smoother;
sigma_smoother = Bayes.SIGMA_smoother;

K = length(mu_smoother);

%% likelihood of state process 
state_val = 0 ;
for k=2 : K
    x_k_k_1 =  mu_smoother(:,k) - mu_smoother(:,k);
  
    state_val = state_val + (x_k_k_1.' * Q * x_k_k_1);
end
likeLihood_val = state_val + (K-1)* log(det(Q));

end

%% nested function 
