function [MU_filter_tot, SIGMA_filter_tot, MU_oneStep_tot, SIGMA_oneStep_tot] = gc_bayes_filter_oneStep_switch(Y_k, L, mu, incongruent_vec, W, mu_1_0, sigma_1_0, F, D, Q, dim, mat_switch_ind)
% GC_BAYES_FILTER this function estimate the posterior probability of the
% state by estimating the Likelihood and One-Step Prediction distributions
% (for further information please take a look at part 1-5 of supplemntary)
%
% INPUTS
% f_xx      : state transient matrix
% p_xx_pre  : distribution of previous posterior
% y_k       : the current observation
% x         : an interval of state variables
% L         : given eigenvector
% mu        : given mu
% incongruent_vec : 1 trail is incngruent, 0 trial is congruent
% W               : matrix of weight
% mu_1_0          : initial value for mu of filtre used in GA
% sigma_1_0       : initial value for sigma of filter used in GA


switch dim
    case 1
        mu_1_0 = mu_1_0(1,1);
        sigma_1_0 = sigma_1_0(1,1);
        F = F(1,1);
        D = D(1,1);
        Q = Q(1,1);
    case 2

end

l_tot = length(Y_k);

MU_filter_tot = zeros(1, l_tot);
SIGMA_filter_tot = zeros(l_tot, 1);
MU_oneStep_tot = zeros(1, l_tot);
SIGMA_oneStep_tot = zeros(l_tot, 1);


vec_active_switch = mat_switch_ind(1, :);

Num_seg = length(vec_active_switch);

seg_ = 0;
while seg_<Num_seg
    seg_ = seg_ + 1;

    active_status = vec_active_switch(seg_);
    
    if active_status==1
        
        ind_l = mat_switch_ind(2, seg_);
        ind_u = mat_switch_ind(3, seg_);

        Y_k_seg = Y_k(ind_l: ind_u, :);

        if seg_>1
            ind_u_seg_1 = mat_switch_ind(3, seg_-1);
            mu_0_seg = MU_filter_tot(1, ind_u_seg_1);
            mu_1_0 = mu_0_seg;
        end


        m = size(Y_k_seg);
        % m(1) = K (number of states)
        % m(2) = number of channels
        K = m(1);



        MU_filter = zeros(dim, K);
        SIGMA_filter = zeros(K, dim, dim);
        MU_oneStep = zeros(dim, K);
        SIGMA_oneStep = zeros(K, dim, dim);

        MU_oneStep(:, 1) = mu_1_0;
        SIGMA_oneStep(1, :, :) = sigma_1_0;


        % % use squeeze(SIGMA_FILTER(1, :, :)) to convert 3d into 2d

        % estimation of posterior probablity of x based BAYES FILTER
        for i=1:K
            %     disp(i)
            y_k = Y_k_seg(i,:).';  % current observation


            mu_k_1 = MU_oneStep(:, i);
            sigma_k_1 = squeeze(SIGMA_oneStep(i, :, :));

            % calculate U
            %     Ly = L*(y_k-mu);
            Ly = conj(L')*(y_k-mu);
            U = Ly.*conj(Ly);


            m = size(Ly);



            %% estimation of sigma
            switch dim
                case 1
                    sigma_W_lambda_U_m = zeros(dim);
                    for j=1 : m(1)
                        w_m = W(j, 1);
                        w_bias = W(j,2);
                        lambda_m = exp(mu_k_1(1,1) * w_m(1,1) + w_bias);
                        sigma_W_lambda_U_m = sigma_W_lambda_U_m + w_m.^2 * inv(lambda_m) * U(j, 1);
                    end


                case 2
                    sigma_W_lambda_U_m = zeros(dim);
                    for j=1 : m(1)
                        w_m = W(j, 1:2);
                        w_bias = W(j,3);
                        lambda_m = exp(mu_k_1(1,1) * w_m(1,1) + incongruent_vec(i, 1) *  mu_k_1(2, 1) * w_m(1, 2) + w_bias);
                        sigma_W_lambda_U_m = sigma_W_lambda_U_m + w_m.' * w_m * inv(lambda_m) * U(j, 1);
                    end
            end


            SIGMA_FILTER_inv = inv(sigma_k_1) + sigma_W_lambda_U_m;
            sigma_filter = inv(SIGMA_FILTER_inv);
            SIGMA_filter(i, :, :) = sigma_filter;

            %% estimation of mu
            switch dim
                %1 d
                case 1
                    sigma_task = zeros(m(1), m(1));
                    for j=1: m(1)
                        w_m = W(j, 1);
                        w_bias = W(j,2);
                        lambda_m = exp(mu_k_1(1,1)*w_m(1,1) + w_bias);
                        sigma_task(j, j) = lambda_m;
                    end
                    W_temp = W(:, 1);


                    % 2d
                case 2
                    sigma_task = zeros(m(1), m(1));
                    for j=1: m(1)
                        w_m = W(j, 1:2);
                        w_bias = W(j,3);
                        lambda_m = exp(mu_k_1(1,1)*w_m(1,1) + incongruent_vec(i, 1)* mu_k_1(2, 1)*w_m(1, 2) + w_bias);
                        sigma_task(j, j) = lambda_m;
                    end
                    W_temp = W(:, 1:2);

            end


            mu_filter = mu_k_1 + sigma_filter * (W_temp.' * inv(sigma_task) * U - W_temp.' * ones(m(1),1));
            MU_filter(:, i) =  mu_filter;


            mu_k_1 = MU_filter(:, i);
            sigma_k_1 = squeeze(SIGMA_filter(i, :, :));

            %% update oneStep
            MU_oneStep(:, i+1) = F*mu_k_1 + D;
            SIGMA_oneStep(i+1, :, :) = F*sigma_k_1*F.' + Q;
        end
        SIGMA_oneStep(i+1, :, :) = [];




    elseif active_status==0

        stop_= 0;
        cnt_ = 1;
        while stop_==0
            seg_i = seg_+cnt_;

            if seg_i<Num_seg
                state_seg_next = vec_active_switch(seg_i);


                diff_ = abs(active_status-state_seg_next);
                if diff_ == 0
                    cnt_ = cnt_ + 1;
                else
                    stop_ = 1;
                end
            else
                stop_ = 1;
            end
        end

        seg_end = seg_ + cnt_ -1;

        ind_l = mat_switch_ind(2, seg_);
        ind_u = mat_switch_ind(3, seg_end);
        
        Y_k_seg = Y_k(ind_l: ind_u, :);



        if seg_>1
            ind_u_seg_1 = mat_switch_ind(3, seg_-1);
            mu_0_seg = MU_filter_tot(1, ind_u_seg_1);
        else
            mu_0_seg =0;
        end


        K_ = length(Y_k_seg);
        MU_x = zeros(1, K_);
        mu_filt = zeros(1, K_);
        sigma_filt = Q*ones(K_, 1);
        sigma_smoother = Q*ones(K_, 1);

        mu_filt(1,1) = mu_0_seg;

        cnt_smoother = 0;
        for kk=2: K_
            MU_x(:, kk) = F*mu_filt(:, kk-1) + D; % mu of normal distribution

            mu_filt(:, kk) = normrnd(MU_x(:, kk) ,Q);
            sigma_filt(kk, 1) = kk*Q;

            if kk< (K_/2 + 5)
                sigma_smoother(kk, 1) = kk*Q;
            else
                cnt_smoother = cnt_smoother +2*normrnd(1, 0.4);
                sigma_smoother(kk, 1) = (kk-cnt_smoother)*Q;
            end
        end


        MU_filter = mu_filt;
        SIGMA_filter = sigma_filt;
        MU_oneStep = MU_filter;
        SIGMA_oneStep = SIGMA_filter;

        seg_ = seg_end;
    end


    MU_filter_tot(1, ind_l: ind_u) = MU_filter;
    SIGMA_filter_tot(ind_l: ind_u, 1) = SIGMA_filter;
    MU_oneStep_tot(1, ind_l: ind_u) = MU_oneStep(1, 1: length(MU_filter));
    SIGMA_oneStep_tot(ind_l: ind_u, 1) = SIGMA_oneStep;
    

end

end



