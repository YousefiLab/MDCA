function [gc_val] = gc_pdf_gc(MU_smoother, SIGMA_smoother, W)
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

mw = size(W);
dim = mw(2)-1;

m = size(MU_smoother);
K=m(2);

gc_val = zeros(3, K);
for k =1: m(2)
    for kk=1:3
        
        switch kk
            case 1 % middle band
                x_temp = MU_smoother(k);
                sigma_k = fun_sigma_calc(x_temp, W, dim);
                diag_sigma = diag(sigma_k);
                gc_val(1, k) = max(diag_sigma)./sum(diag_sigma);
                
            case 2 %upper band
                x_temp = MU_smoother(k)+ 2*sqrt(SIGMA_smoother(k));
                sigma_k = fun_sigma_calc(x_temp, W, dim);
                diag_sigma = diag(sigma_k);
                gc_val(2, k) = max(diag_sigma)./sum(diag_sigma);
                
            case 3 % lower band
                x_temp = MU_smoother(k)- 2*sqrt(SIGMA_smoother(k));
                sigma_k = fun_sigma_calc(x_temp, W, dim);
                diag_sigma = diag(sigma_k);
                gc_val(3, k) = max(diag_sigma)./sum(diag_sigma);   
        end
    
    
    end
end




end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global Coherence estimator
function [Sigma_W_x] = fun_sigma_calc(x_temp, W, dim)
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

mw = size(W);
ch_num  = mw(1);

switch dim
    
    %1 d
    case 1
        Sigma_W_x = zeros(ch_num);
        for j=1: ch_num
            w_m = W(j, 1);
            w_bias = W(j,2);
            lambda_m = exp(x_temp(1,1) * w_m(1,1)+ w_bias);
            Sigma_W_x(j, j) = lambda_m;
        end
        
        
        % 2d
    case 2
        Sigma_W_x = zeros(ch_num);
        for j=1: ch_num
            w_m = W(j, 1:2);
            w_bias = W(j,3);
            lambda_m = exp(x_temp(1,1) * w_m(1,1)  +   incongruent_vec(i, 1) *  x_temp(2, 1) * w_m(1, 2) + w_bias);
            Sigma_W_x(j, j) = lambda_m;
        end     
end

end



