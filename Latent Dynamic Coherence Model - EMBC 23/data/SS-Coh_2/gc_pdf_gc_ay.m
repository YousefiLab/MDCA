function [pmf_gc_tot, s] = gc_pdf_gc(x, mu, param, dist_prob, step_gc)
% GC_PDF_X calculates pdf of global coherence for a given distribution
% probability
%
% INPUTS
% x_k        : state variables
% x          : interval of x used for estimated dist_prob
% g_k        : estimated Global Coherence by the proposed model
% dist_prob  : probability distribution for all x and different steps
% step_gc    : the step of is used to calculate GC's pdf
%
% OUTPUTS
% pdf_gc     : pdf of Global Coherence

if nargin < 4
    error('4 first inputs "x, mu, param, dist_prob" should be passed')
end
if nargin== 4
    step_gc = .001;
end

gc_x = gc_gk_estimator(x, mu, param);

[gc_x_sorted, I] = sort(gc_x);

% interval of estimated GC
s = 0:step_gc:1;

m = size(dist_prob);
% m(1) = length(x)
% m(2) = K

% pmf of all steps
pmf_gc_tot = zeros(length(s), m(2));

for k =1: m(2)
    prob_xk = dist_prob(:,k);
    
    % CDF probability of GC for step k
    cdf_gc_k = prob_xk(I);
    res      = zeros(length(s),1);
    for t=1:length(I)
        [~,t_ind]  = min(abs(gc_x_sorted(t)-s));
        if t_ind > 1 || t_ind <length(s) 
            res(t_ind) = res(t_ind)+cdf_gc_k(t);    
        else
            res(t_ind-1) = res(t_ind)+0.25*cdf_gc_k(t);
            res(t_ind) = res(t_ind)+0.5*cdf_gc_k(t);
            res(t_ind+1) = res(t_ind)+0.25*cdf_gc_k(t);
        end
    end
    
    % gc_x_sorted;
    
    
    
    % CDF probability of GC for step k
    
    
%     for i= 1:length(s)
%         temp_cdf = 0;
%         % indices that we want to evaluate GC value on them (shows in which step (k) we have desired GC)
%         
%         ind_before = find(gc_x_sorted<=s(i));
%         
%         if ~isempty(ind_before)
%             temp_cdf = interpolate_fun(ind_before, gc_x_sorted, s(i), I, prob_xk);
%         end
%         
%         cdf_gc_k(1,i) = temp_cdf;
%     end
%     
%     % normalizing CDF
%     cdf_gc_k = cdf_gc_k./cdf_gc_k(end);
%     
%     % getting pdf from CDF for step k
%     pmf_gc_k = zeros(1,length(s));
%     
%     for i=2:length(s)
%         pmf_gc_k(1,i) = cdf_gc_k(1,i) - cdf_gc_k(1,i-1);
%     end
%     
    pmf_gc_tot(:,k) = res;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global Coherence estimator
function [g_k] = gc_gk_estimator(x, mu, param)
% GC_GK_ESTIMATOR estimate Global Coherence for given interval of x
%
% INPUTS
% x        : interval of x that we want to evaluate
% mu       : added mean to estimated observatoin
% param.a  : a value for generating lambda (used in eigenvalue matrix D)
% param.b  : b value for generating lambda (used in eigenvalue matrix D)
%
% OUTPUTS
% g_k      : estimate the level of coherence at each time

if nargin < 3
    error('all inputs "x_k" should be passed')
end

n_p = length(x);    % Number of Points
d_p = length(mu);   % Dimension of the Problem


% lambda_k : contains the values of lambda that is used at each step for generating D_k
% eigenvalues that are generated in each time step (k)
lambda_k = zeros(n_p , d_p);

% this vecotor contains the GC value based proposed model approach
g_k = zeros(n_p,1);

for k=1 : n_p
    
    % calculate D for x_k and given am, bm
    D_k = diag(exp(param.a + x(k).*param.b));  % D
    
    % store values of lambda at each time
    lambda_k(k,:) = diag(D_k);
    
    % calculate g_k
    g_k(k,1) = (max(lambda_k(k,:).')./sum(lambda_k(k,:).')).';
end


end

%% interpolate
function [temp_cdf] = interpolate_fun(ind_before, gc_x_sorted, s_index, I, prob_xk)


ind_before_correspond = I(ind_before);

ind_after = min(find(gc_x_sorted> s_index));
if  isempty(ind_after)
    
    temp_cdf = sum(prob_xk(ind_before_correspond));
    temp_cdf2 = sum(prob_xk(ind_before_correspond));
else
    ind_after_correspond = I(ind_after);
    ind_bef_aft = [ind_before_correspond; ind_after_correspond];
    
    gc_bef = gc_x_sorted(ind_before(end));
    gc_aft = gc_x_sorted(ind_after);
    
    
    temp_cdf_before = sum(prob_xk(ind_before_correspond)); % calculate cdf for specified GC value
    temp_cdf_after = sum(prob_xk(ind_bef_aft));
    
    
    x = [gc_bef(end) gc_aft];
    y = [temp_cdf_before temp_cdf_after];
    
    % interpolate first order
%     x
%     y
    temp_cdf1 = interp1(x,y,s_index);
    
    % interpolate second order
    temp_cdf2 = spline(x, y, s_index);
end

temp_cdf = temp_cdf2;
end


