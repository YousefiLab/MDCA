function [MU_smoother_tot, SIGMA_smoother_tot, G_k_tot]= gc_bayes_smoother_switch(MU_filter_tot, SIGMA_filter_tot, MU_oneStep_tot, SIGMA_oneStep_tot, F, Q, dim, mat_switch_ind)
% GC_BAYES_SMOOTHER used to estimate the posterior probability based Bayes

switch dim
    case 1
        F = F(1,1);

    case 2

end

l_tot = length(MU_filter_tot);

MU_smoother_tot = zeros(1, l_tot);
SIGMA_smoother_tot = zeros(l_tot, 1);
m_sigma = size(SIGMA_filter_tot);
G_k_tot = zeros(l_tot, 1);


vec_active_switch = mat_switch_ind(1, :);

Num_seg = length(vec_active_switch);

seg_ = 0;
while seg_< Num_seg
    seg_ = seg_ + 1;
    active_status = vec_active_switch(seg_);
    
    if active_status==1


        ind_l = mat_switch_ind(2, seg_);
        ind_u = mat_switch_ind(3, seg_);


        MU_filter = MU_filter_tot(1, ind_l: ind_u);
        SIGMA_filter = SIGMA_filter_tot(ind_l: ind_u, 1);
        MU_oneStep = MU_oneStep_tot(1, ind_l: ind_u);
        SIGMA_oneStep = SIGMA_oneStep_tot(ind_l: ind_u, 1);




        m_mu = size(MU_filter);
        m_sigma = size(SIGMA_filter);
        K = m_mu(2);

        MU_smoother = zeros(m_mu(1), m_mu(2));
        MU_smoother(:, K)= MU_filter(:, K);

        switch dim
            % d1
            case 1
                SIGMA_smoother = zeros(m_sigma(1), m_sigma(2));
                G_k = zeros(m_sigma(1), m_sigma(2));
                SIGMA_smoother(K, 1)= SIGMA_filter(K, 1);

                % d2
            case 2
                SIGMA_smoother = zeros(m_sigma(1), m_sigma(2), m_sigma(3));
                G_k = zeros(m_sigma(1), m_sigma(2), m_sigma(3));

                SIGMA_smoother(K, :, :)= squeeze(SIGMA_filter(K, :, :));

        end



        for i=K-1:-1: 1
            mu_filter_curr = MU_filter(:, i);
            sigma_filter_curr = squeeze(SIGMA_filter(i, :, :));

            mu_oneStep_next = MU_oneStep(:, i+1);
            sigma_oneStep_next = squeeze(SIGMA_oneStep(i+1, :, :));

            mu_smoother_next = MU_smoother(:, i+1);
            sigma_smoother_next= squeeze(SIGMA_smoother(i+1, :, :));

            Gk = sigma_filter_curr*F.'*inv(sigma_oneStep_next);
            G_k(i, :, :) = Gk;

            mu_smoother_curr = mu_filter_curr + Gk*(mu_smoother_next - mu_oneStep_next);
            simga_smoother_curr = sigma_filter_curr + Gk* (sigma_smoother_next- sigma_oneStep_next)*Gk';

            MU_smoother(:, i)= mu_smoother_curr;
            SIGMA_smoother(i, :, :)= simga_smoother_curr;


        end

        

        MU_smoother_tot(1, ind_l:ind_u) = MU_smoother;
        SIGMA_smoother_tot(ind_l:ind_u, 1) = SIGMA_smoother;

        G_k(end, :) = G_k(end-1, :);
        G_k_tot(ind_l:ind_u, 1) = G_k;

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

        MU_filter = MU_filter_tot(1, ind_l: ind_u);
        SIGMA_filter = SIGMA_filter_tot(ind_l: ind_u, 1);
        MU_oneStep = MU_oneStep_tot(1, ind_l: ind_u);
        SIGMA_oneStep = SIGMA_oneStep_tot(ind_l: ind_u, 1);
        Gk_ = zeros(length(MU_filter), 1);

        K_ = length(MU_filter);
        sigma_smoother = Q(1,1)*ones(K_, 1);


        cnt_smoother = 0;
        for kk=2: K_
%             str_ = sprintf('%d', kk);
%             disp(str_)
            sigma_filter_curr = SIGMA_filter(kk);
            sigma_oneStep_next = SIGMA_oneStep(kk);
            Gk_temp = sigma_filter_curr*F.'*inv(sigma_oneStep_next-1);
            if isnan(Gk_temp)
                Gk_temp = 0; 
            end
            Gk_(kk-1, 1) = Gk_temp;
            
            if kk< (K_/2 + 5)
                sigma_smoother(kk, 1) = kk*Q;
            else
                cnt_smoother = cnt_smoother +2*normrnd(1, 0.4);
                sigma_smoother(kk, 1) = (kk-cnt_smoother)*Q;
            end
        end

        MU_smoother_tot(1, ind_l:ind_u) = MU_filter;
        SIGMA_smoother_tot(ind_l:ind_u, 1) = sigma_smoother;

        Gk_(end, :) = Gk_(end-1, :);
        G_k_tot(ind_l:ind_u, 1) = Gk_;
        
%         figure()
%         plot(sigma_smoother)

        seg_ = seg_end;
    end
end

end


