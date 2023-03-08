function [f_xx] = gc_state_transition(x,s_v,a,b)
% GC_STATE_TRANSIENT calculates state transient matrix for different values
% of x based a normal distribution - p(x(k)|x(k-1))
%
% INPUTS
% x      : an interval of all possible state variables
% s_v    : given var
% a      : a paramater for generating mu
% b      : b paramteter for genrating mu
%
% OUTPUTS
% f_xx   : state transient matrix

if nargin == 0
    error('the first input "x" should be passed')
end
if nargin == 1
    s_v = .1;
    a = 1;
    b = 0;
end
if nargin == 2
    a = 1;
end
if nargin == 3
    b = 0;
end

% transient state matrix
f_xx = zeros(length(x) , length(x));

for i=1:length(x)
    
    m_0 = a*x(i) + b;         % mu of normal distribution
    std_0 = sqrt(s_v);        % variance of normal distribution
    
    % calculate normal pdf for elemants of x with a given mu (x_0) and var (s_v)
    f_xx(:,i) = normpdf(x,m_0,std_0);
    f_xx(:,i) = f_xx(:,i)/sum(f_xx(:,i));
end

end