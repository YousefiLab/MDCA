function [x_k,x_0] = gc_random_walk(K,s_v,m_0,s_0,a,b)
% GC_RANDOM_WALK generate a sequence of state variables 
% random walk model is used for generating state variables 
% (for further information please take a look at part 1-1 of supplemntary)
%
% INPUTS
% K     : number of states
% s_v   : variance of additive noise in AR model
% m_0   : distribution's mean for generating x(0)
% s_0   : distribution's variance for generating x(0)
% a     : coeficient of x(k) in AR model
% b     : bias value in AR model
%
% OUTPUTS
% x_k   : a sequence of state variables 
% x_0   : initial value used in AR model for generating x_k 


% Initializing (assign default values)
% a: AR(1)      (default 1)
% b: bias       (default 0)
if nargin==0
    a = 1;
    b = 0;
    s_v = 0.01;
    K = 100;
    m_0 = 0;
    s_0 = 0.01;
end
if nargin==1
    a = 1;
    b = 0;
    s_v = 0.01;
    m_0 = 0;
    s_0 = 0.01;
end
if nargin==2
    a = 1;
    b = 0;
    m_0 = 0;
    s_0 = 0.01;
end
if nargin==3
    a = 1;
    b = 0;
    s_0 = 0.01;
end
if nargin==4
    a = 1;
    b = 0;
end
if nargin==5
    b = 0;
end


% Genrating x_0 with a normal distribution with given mean and variance
x_0 = m_0 + randn() * sqrt(s_0);

% x_k is the state variables
x_k = zeros(K,1);

% Generating x_k with the AR model
x_k(1) = a * x_0 + b + randn() * sqrt(s_v);

for k = 2: K
    x_k(k) = a * x_k(k-1) + b + randn() * sqrt(s_v);   % eq (1)
end

end


