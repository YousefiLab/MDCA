function [Y_k,g_k] = gc_observation_simulator(x_k,mu,L,param)
% GC_OBSERVATION_SIMULATOR generate a matrix of samples based on proposed
% model
% Also this function calculates Global Coherence measure
% of model based lambda values of different times
% (for further information please take a look at part 1-2, 1-3 and 1-4 of supplemntary)
%
% INPUTS
% x_k      : random walk variable
% mu       : added mean to estimated observatoin
% L        : Eigenvector
% param.a  : a value for generating lambda (used in eigenvalue matrix D)
% param.b  : b value for generating lambda (used in eigenvalue matrix D)
% 
% OUTPUTS
% Y_k      : observation samples
% g_k      : estimate the level of coherence at each time

if nargin==0
    error('the first input "x_k" should be passed')
end

n_p = length(x_k);  % Number of Points 
d_p = length(mu);   % Dimension of the Problem

% Initializing (assign default values)
if nargin==2
    L = randn(d_p,d_p);
    param.a = .1*randn(length(mu),1);
    param.b = .1*randn(length(mu),1);
end
if nargin==3
    param.a = .1*randn(length(mu),1);
    param.b = .1*randn(length(mu),1);
end

% Y_k is the out put and 
% number of rows is equal to number of points
% number of columns is equal to number of channels
Y_k = zeros(n_p , d_p);

% lambda_k : contains the values of lambda that is used at each step for generating D_k
% eigenvalues that are generated in each time step (k)
lambda_k = zeros(n_p , d_p);

% this vecotor contains the GC value based proposed model approach
g_k = zeros(n_p,1);

Q_k = cell(1, n_p);
for k=1 : n_p
    
    % calculate D for x_k and given am, bm
    D_k = diag(exp(param.a + x_k(k).*param.b));  % D
    D_k_2 = sqrt(D_k);                           % D^.5
 
    % store values of lambda at each time
    lambda_k(k,:) = diag(D_k);                                % eq(7)
    
    % generate a multivariate normal with mu=0 and cov = I
    n_k = mvnrnd(zeros(d_p,1) , eye(d_p,d_p)).';              % eq (4)
    
    % imposing structure of covariance matrix (Q_k(x_k))
    z_k = L*D_k_2*n_k;                                        % eq (5)
    y_k = mu + z_k;                                           % eq (6)
    
    Y_k(k,:) = y_k.';
    
    % calculate g_k 
    g_k(k,1) = (max(lambda_k(k,:).')./sum(lambda_k(k,:).')).'; % eq (8)
end


end