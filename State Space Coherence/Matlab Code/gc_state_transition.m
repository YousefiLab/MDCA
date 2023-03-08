function [f_xx] = gc_state_transition(x, F, D, Q)
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

%     mm_0 = F*x(:,i)+Q;
   
m = size(x);
% transient state matrix
f_xx = zeros(m(1) , m(1));
for i=1: m(1)
    

    
    MU_0 = F*x(i, :)' + D; % mu of normal distribution
    
    % calculate normal pdf for elemants of x with a given mu (x_0) and var (s_v)
    TEMP = mvnpdf(x, MU_0.', Q);    
    f_xx(:,i) = TEMP/sum(TEMP);
    
    %     figure,
%     scatter3(x(:,1),x(:,2),TEMP)

    
%     figure, plot(f_xx(:,i))
end

end