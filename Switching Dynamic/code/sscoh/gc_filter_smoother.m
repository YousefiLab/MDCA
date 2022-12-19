function [Bayes] = gc_filter_smoother(Y_k, Param)
% GC_FILTER_SMOOTHER estimates posterior with bayes filter and smoother ,
% after calculating smoother, we estimate a te distribution for GC
%
% INPUTS
% Y_k                 : observation (measured of FFT)
% x_k                 : state variables
% Param               : is a structure that contains parameters
%    Param.sv         : variance 
%    Param.s_ab       : a and b of random walk model
%    Param.mu         : mean
%    Param.L          : eigenvalue
%    Param.ab         : a structure of am and bm
%          Param.ab.a : am
%          Param.ab.b : bm
% 
% OUTPUTS
% Bayes                    : a structure of bayes filter-smoother results
%     Bayes.smoother       : smoother result 
%     Bayes.filter         : filter result 
%     Bayes.oneStep        : one step prediction of fiter
%     Bayes.filter_prior   : prior probability of filter 
%     Bayes.stateTransient : state transient (step k to k+1)t 
%     Bayes.x              : an interval of x that we want to evalueate
%     Bayes.Y_k            : observation (measured of FFT)
% pmf_gc                   : distribution of Global Coherence
% s                        : an interval for calculating GC


% state transient probabilities - p(x(k)|x(k-1))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % a = Param.s_ab(1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % b = Param.s_ab(2);
dim = Param.dim;
mu = Param.mu;
L = Param.L ;
incongruent_vec = Param.incongruent_vec;
W = Param.W;
mu_1_0 = Param.mu_1_0;
sigma_1_0 = Param.sigma_1_0;
F = Param.F;
D = Param.D;
Q = Param.Q;



% Bayes filter/one step
[MU_filter, SIGMA_filter, MU_oneStep, SIGMA_oneStep] = gc_bayes_filter_oneStep(Y_k, L, mu, incongruent_vec, W, mu_1_0, sigma_1_0, F, D, Q, dim); % 15.a, 15.b


% Bayes smoother
[MU_smoother, SIGMA_smoother, G_k]= gc_bayes_smoother(MU_filter, SIGMA_filter, MU_oneStep, SIGMA_oneStep, F, Q, dim);

% pmf distribution of Global Coherence
GC_val = gc_pdf_gc(MU_smoother, SIGMA_smoother, W);

% structure of bayes filter-smoother results
Bayes.MU_filter = MU_filter;
Bayes.SIGMA_filter = SIGMA_filter;
Bayes.MU_oneStep = MU_oneStep;
Bayes.SIGMA_oneStep = SIGMA_oneStep;
Bayes.MU_smoother = MU_smoother;
Bayes.SIGMA_smoother = SIGMA_smoother;
Bayes.G_k = G_k;
Bayes.GC_val = GC_val;

end