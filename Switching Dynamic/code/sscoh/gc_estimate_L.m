function [L_update, L_diff] = gc_estimate_L(MU_smoother, SIGMA_smoother, Yk, mu, W, L, incongruent_vec, dim)
% GC_L_ESTIMATE estimates the eigenvector matrix (L)


L_init = L;

m = size(L_init);
M = m(1);

L_estimate = zeros(M);

for m=1: M
    
    Wm = W(m, :);
    Wm = Wm.';
    
    % estimate Zm - Eq 7 (UPDATE_L_W.docx)
    Zm = nested_Zm(MU_smoother, SIGMA_smoother, Yk, mu, Wm, dim);
    
    % find Eigenvlaue, Eigenvector
    [eVec,eVal] = eig(Zm);
    
    eVal = diag(abs(eVal));
    
    ind_null = find(eVal==0);
    
    if isempty(ind_null)
        ind_min_eVal = find(eVal == min(eVal));
        eVec_ans = eVec(:, ind_min_eVal);
    else
        eVec_ans = eVec(:, ind_null);
    end
    
    Um_estimate = eVec_ans;
    L_estimate(:, m) = Um_estimate;
    
    Um_init = L_init(:, m);
end

L_estimate_orth = orth(L_estimate);

L_update = L_estimate_orth;

L_diff = L_update - L_init;

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Zm] = nested_Zm(MU_smoother, SIGMA_smoother, Yk, mu, Wm, dim)

K = length(MU_smoother);
M = length(mu);

switch dim
    case 1
        Wm_hat = Wm(1,1);
    case 2
        Wm_hat = Wm(1:2);
end
Wm_biased = Wm(3, 1);

Zm = zeros(M);

for k=1: K
    % p_k - Eq 2.a (UPDATE_L_W.docx)
    y_k = Yk(k, :).';
    p_k = y_k- mu;
    
    mu_s = MU_smoother(:, k);
    sigma_s = squeeze(SIGMA_smoother(k, :, :));
    
    % estimate E[lambda_inv] - Eq 4.c (UPDATE_L_W.docx)
    E_lambda_inv = exp((-Wm_hat.' * (mu_s-.5*sigma_s*Wm_hat)) - Wm_biased);

    Zm = Zm + E_lambda_inv*p_k*p_k';
end

end



