clc
clear
close all

dim = 1; % 1 or 2
ch_num = 5; % 5, 8, 10, 20
Iter = 11;
num_switch = 1;


% plt_ml = 0;
% if plt_ml==1
%     ml = [243, 350, 500, 600, 760, 800, 1100. 1215, 1502, 1590, 1592, 1592.5, 1593, 1593.5, 1600, 1600];
%     figure()
%     plot(ml,'k', 'linewidth', 4)
%     xlabel('Iteration')
%     ylabel('ML Growth (max Q)')
%     ylim([min(ml)-100, max(ml)+100])
%     xlim([1, length(ml)])
%     ax = gca;
%     ax.FontSize = 16;
% end

use_bias = 0; % if 1, use all bias elements  for W
if use_bias==0
    use_bias_2_end = 1; % if 1 use all bias elements except first
end

using_simulation = 0; % 1: using sim, 0: using real data
if using_simulation==0

    coef_noise_W = 0.5;
    coef_noise_mu = 0;

    coef_w = 1;

    coef_length_sim = 1;

    new_val = 1;
    new_W = new_val;
    new_SIM_data = new_val;

    fprintf('Using Real Data Set for Analysing SS-GCoh Method')
    fprintf('\n------------------------------\n')

elseif using_simulation==1
    coef_noise_W = 1;

    coef_noise_mu = 3;

    coef_w = .1;

    coef_length_sim = .1;

    new_val = 1;
    new_W = new_val;
    new_SIM_data = new_val;

    fprintf('Using Simulation Data Set for Analysing SS-GCoh Method')
    fprintf('\n------------------------------\n')
end


col_W_1 = 2: ch_num;
col_W_2 = 2:ch_num;
if use_bias ==1
    col_W_3 = 1: ch_num;

elseif use_bias_2_end==1
    col_W_3 = 2: ch_num;
end

QQ_sim = [15 0; 0 80].*10^-3;

%% PARAMETER UPDATE
ParamaUpdate.sigma_update = 0; % 1 update, 0 don't update
ParamaUpdate.mu_update = 0;
ParamaUpdate.L_update  = 0;
ParamaUpdate.W_update  = 1;

% based on switch we need to run the model K times and then connect them to
% each other

if num_switch==3
    structure_switch = [1, 2, 3, 2, 1];
elseif num_switch==2
    structure_switch = [1, 2, 1];
end

Num_seg = length(structure_switch);

PARAM_q1= cell(Num_seg, Iter+1);
PARAM_q2= cell(Num_seg, Iter+1);
Bayes_Step_q1 = cell(Num_seg, Iter);
Bayes_Step_q2 = cell(Num_seg, Iter);
if num_switch==3
    PARAM_q3= cell(Num_seg, Iter+1);
    Bayes_Step_q3 = cell(Num_seg, Iter);
end

for q_k=1: num_switch-1


    % load inilial points
    [cfg_Init, Y_k]  = gc_load_yk_cfg(ch_num, num_switch);
    
    l_tot = length(Y_k);
    if q_k ==1
        gc_print_fun_1(cfg_Init, ch_num) % print some info
    end

    

    %% calc switch information
    empirical_switch_point = cfg_Init.empirical_switch_pionts;
    vec_1 = [empirical_switch_point, 0];
    vec_2 = [0, empirical_switch_point];

    vec_diff = abs(vec_1-vec_2);
%     figure()
%     stem(vec_diff)

    ind_jump = find(vec_diff>0);
    ind_jump(end) = ind_jump(end);
    mat_switch_ind = zeros(3, length(structure_switch));
    mat_switch_ind(1, :) = structure_switch;

    for ii=1: length(structure_switch)
        mat_switch_ind(2, ii) = ind_jump(ii);
        mat_switch_ind(3, ii) = ind_jump(ii+1)-1;
    end


    %% end of index
    

    %% segment each state and run one state model for it 


    % generate new W
     if new_W==1
        W = gc_gen_W_init(coef_w, ch_num);
    end

    L_init = cfg_Init.L;
    L_init_switch = L_init{q_k};
    

    
    for seg_=1: Num_seg
        % print state model
        str_print = sprintf('State Space Model = %d/%d  - seg = %d/%d', q_k, num_switch, seg_, Num_seg);
        disp(str_print)


        switch_active = 0;
        if structure_switch(seg_) == q_k
            switch_active =1;
        end


        ind_strt_seg = mat_switch_ind(2, seg_);
        ind_end_seg =  mat_switch_ind(3, seg_);

        
        Yk_seg = Y_k(ind_strt_seg: ind_end_seg, :);
        t_min = cfg_Init.cfgGen.t_min;



        
        y_end = length(Yk_seg);

       

        K = length(Yk_seg);         % number of all stepes
        KK = floor(coef_length_sim*K);
        KK = 1000;

        if new_SIM_data == 1
            [xk_sample_tot, Yk_sample_tot, cnfg_sim_tot] = gc_sim_data(dim, QQ_sim, L_init, KK, Y_k, ch_num, q_k);
            [xk_sample, Yk_sample, cnfg_sim] = gc_sim_data(dim, QQ_sim, L_init, KK, Yk_seg, ch_num, q_k);
            %         figure, plot(xk_sample)
        else
            str_load = sprintf('X_Y_sim%d_good.mat', dim);
            load(str_load)
            %     figure, plot(xk_sample)
        end


        %% NEW INITIAL PARAMS
        [init_Param, update_coef_W, W_sim] = gc_init_param(ParamaUpdate, cnfg_sim, coef_noise_mu, col_W_1, col_W_2, col_W_3, coef_noise_W, ch_num, dim);

        W_init = init_Param.W;

        %% init param
        init_Param.incongruent_vec = cnfg_sim.incongruent_vec;
        init_Param.mu_1_0 = [0; 0];
        init_Param.sigma_1_0 = 300*[1 0; 0 1]*10^-3;
        init_Param.F = cnfg_sim.F;
        init_Param.D = cnfg_sim.D;
        init_Param.dim = cnfg_sim.dim;
        init_Param.L = L_init_switch;
        
        %% initialize Param and Bayes container 
        if num_switch ==2

            if q_k==1
                PARAM_q1{seg_, 1}= init_Param;
            elseif q_k==2
                PARAM_q2{seg_, 1}= init_Param;
            end

        elseif num_switch==3
            if q_k==1
                PARAM_q1{seg_, 1}= init_Param;
            elseif q_k==2
                PARAM_q2{seg_, 1}= init_Param;
            elseif q_k==3
                PARAM_q3{seg_, 1}= init_Param;
            end
        end


        %% training model
        
        sv_update = [];

        step_interval = 1:length(Yk_seg);
        break_con = 0;
        thr_q = 10^-5;


        % ParamaUpdate.W_update  = 0;

        if using_simulation == 1
            Yk = Yk_sample;
        end

        % model training section 
        for iter=1:Iter

            disp(iter)
            disp('-----------')
            
            if q_k==1
                param_curr = PARAM_q1{seg_, iter};
            elseif q_k==2
                param_curr = PARAM_q2{seg_, iter};
            end

            


            if switch_active==1
                %%% filter - smoother
                Bayes = gc_filter_smoother(Yk_seg, param_curr);

                 %%% update parameters
                [updated_param, error_tot] = gc_parameter_update(param_curr, ParamaUpdate, Bayes, Yk_seg, iter, W_sim, update_coef_W);
            else
                %%% follow state-transient
                F = 1;
                D = 0;
                Q = 0.015*10e-2;

                K_ = length(Yk_seg);
                MU_x = zeros(1, K_);
                mu_filt = zeros(1, K_);
                sigma_filt = Q*ones(K_, 1);
                sigma_smoother = Q*ones(K_, 1);
                
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
                MU_smoother = MU_filter;
                SIGMA_smoother = sigma_smoother;
                

                

                plt_ = 0; 
                if plt_ ==1
                    figure()
                    subplot(211)
                    plot(MU_filter);
                    hold on
                    lower_up = MU_filter + SIGMA_filter.';
                    lower_dwn = MU_filter - SIGMA_filter.';
                    plot(lower_up)
                    hold on
                    plot(lower_dwn)
                    title('Filter')

                    subplot(212)
                    plot(MU_smoother);
                    hold on
                    lower_up = MU_smoother + SIGMA_smoother.';
                    lower_dwn = MU_smoother - SIGMA_smoother.';
                    plot(lower_up)
                    hold on
                    plot(lower_dwn)
                    title('Smoother')
                end



                W = param_curr.W;
                GC_val = gc_pdf_gc(MU_smoother, SIGMA_smoother, W);

                Bayes= [];
                Bayes.MU_smoother = MU_smoother;
                Bayes.SIGMA_smoother = SIGMA_smoother;
                Bayes.MU_filter = MU_filter;
                Bayes.SIGMA_filter = SIGMA_filter;
                Bayes.MU_oneStep = MU_oneStep;
                Bayes.SIGMA_oneStep = SIGMA_oneStep;
                Bayes.GC_val = GC_val;
                
                % we don't update parameters 
                updated_param = param_curr;
            end


           


%             MU_iter(iter, :) = param_curr.mu;
% 
%             W_iter(iter, :, :) = param_curr.W;
% 
%             % absolute error
%             abs_error_TOT_iter(1, iter) = error_tot;
% 
            
            if q_k==1
                Bayes_Step_q1{seg_, iter} = Bayes;
                PARAM_q1{seg_, iter+1} = updated_param;
            elseif q_k==2
                Bayes_Step_q2{seg_, iter} = Bayes;
                PARAM_q2{seg_, iter+1} = updated_param;
            end

        end
    end
end


% MU_iter(iter+1, :) = updated_param.mu;
% W_iter(iter+1, :, :) = updated_param.W;
% 
% iter_first = 1;
% iter_last =  length(Bayes_Step);
% 
% bayes_first = Bayes_Step{iter_first};
% bayes_last = Bayes_Step{iter_last};
% bayes_last = Bayes_Step{5};
%% here I need to find a way to concatenate three different switches to each other and then pass them for the next step

iter_first = 1;
iter_last = Iter;

% num_switch
for q_k=1: num_switch


    if q_k==1
        Bayes = Bayes_Step_q1;
        PARAM = PARAM_q1;
    elseif q_k==2
        Bayes = Bayes_Step_q2;
        PARAM = PARAM_q2;
    end
    
    
    MU_smoother_concat_first = zeros(1, l_tot);
    SIGMA_smoother_concat_first = zeros(l_tot, 1);
    MU_filter_concat_first = zeros(1, l_tot);
    SIGMA_filterr_concat_first = zeros(l_tot, 1);
    MU_oneStep_concat_first = zeros(1, l_tot);
    SIGMA_oneStep_concat_first = zeros(l_tot, 1);
    GC_val_concat_first  = zeros(3, l_tot);


    MU_smoother_concat_last = zeros(1, l_tot);
    SIGMA_smoother_concat_last = zeros(l_tot, 1);
    MU_filter_concat_last = zeros(1, l_tot);
    SIGMA_filterr_concat_last = zeros(l_tot, 1);
    MU_oneStep_concat_last = zeros(1, l_tot);
    SIGMA_oneStep_concat_last = zeros(l_tot, 1);
    GC_val_concat_last  = zeros(3, l_tot);
    
    ind_strt = 1; 
    for seg_=1: Num_seg
        Bayes_seg_first = Bayes{seg_, 1};
        Bayes_seg_last = Bayes{seg_, iter_last};

        
        
        MU_smoother_seg_first  = Bayes_seg_first.MU_smoother;
        SIGMA_smoother_seg_first  = Bayes_seg_first.SIGMA_smoother;
        MU_filter_seg_first  = Bayes_seg_first.MU_filter;
        SIGMA_filter_seg_first  = Bayes_seg_first.SIGMA_filter;
        MU_oneStep_seg_first  = Bayes_seg_first.MU_oneStep;
        SIGMA_oneStep_seg_first  = Bayes_seg_first.SIGMA_oneStep;
        GC_val_seg_first  = Bayes_seg_first.GC_val;

        MU_smoother_seg_last  = Bayes_seg_last.MU_smoother;
        SIGMA_smoother_seg_last   = Bayes_seg_last.SIGMA_smoother;
        MU_filter_seg_last   = Bayes_seg_last.MU_filter;
        SIGMA_filter_seg_last   = Bayes_seg_last.SIGMA_filter;
        MU_oneStep_seg_last   = Bayes_seg_last.MU_oneStep;
        SIGMA_oneStep_seg_last   = Bayes_seg_last.SIGMA_oneStep;
        GC_val_seg_last   = Bayes_seg_last.GC_val;

        l_seg = length(MU_smoother_seg_first);
        
        ind_end = ind_strt+l_seg-1;

        MU_smoother_concat_first(1, ind_strt: ind_end) = MU_smoother_seg_first;
        SIGMA_smoother_concat_first(ind_strt: ind_end, 1) = SIGMA_smoother_seg_first;
        MU_filter_concat_first(1, ind_strt: ind_end) = MU_filter_seg_first;
        SIGMA_filterr_concat_first(ind_strt: ind_end, 1) = SIGMA_filter_seg_first;
        MU_oneStep_concat_first(1, ind_strt: ind_end) = MU_oneStep_seg_first(1, 1: l_seg);
        SIGMA_oneStep_concat_first(ind_strt: ind_end, 1) = SIGMA_oneStep_seg_first;
        GC_val_concat_first(:, ind_strt: ind_end)  = GC_val_seg_first;
        



        MU_smoother_concat_last(1, ind_strt: ind_end) = MU_smoother_seg_last;
        SIGMA_smoother_concat_last(ind_strt: ind_end, 1) = SIGMA_smoother_seg_last;
        MU_filter_concat_last(1, ind_strt: ind_end) = MU_filter_seg_last;
        SIGMA_filterr_concat_last(ind_strt: ind_end, 1) = SIGMA_filter_seg_last;
        MU_oneStep_concat_last(1, ind_strt: ind_end) = MU_oneStep_seg_last(1, 1: l_seg);
        SIGMA_oneStep_concat_last(ind_strt: ind_end, 1) = SIGMA_oneStep_seg_last;
        GC_val_concat_last(:, ind_strt: ind_end)  = GC_val_seg_last;
        
        ind_strt = ind_end+1;

    end



    bayes_first = [];
    bayes_first.MU_smoother= MU_smoother_concat_first;
    bayes_first.SIGMA_smoother= SIGMA_smoother_concat_first;
    bayes_first.MU_filter= MU_filter_concat_first;
    bayes_first.SIGMA_filter= SIGMA_filterr_concat_first;
    bayes_first.MU_oneStep= MU_oneStep_concat_first;
    bayes_first.SIGMA_oneStep= SIGMA_oneStep_concat_first;
    bayes_first.GC_val= GC_val_concat_first;


    bayes_last = [];
    bayes_last.MU_smoother= MU_smoother_concat_last;
    bayes_last.SIGMA_smoother= SIGMA_smoother_concat_last;
    bayes_last.MU_filter= MU_filter_concat_last;
    bayes_last.SIGMA_filter= SIGMA_filterr_concat_last;
    bayes_last.MU_oneStep= MU_oneStep_concat_last;
    bayes_last.SIGMA_oneStep= SIGMA_oneStep_concat_last;
    bayes_last.GC_val= GC_val_concat_last;


    % finding the y range
    strct_min_max = fun_plt_find_yrange(bayes_first, bayes_last, xk_sample_tot, dim, using_simulation);

    %% plot sim
    % fun_plt_simulated_data(xk_sample, dim, strct_min_max)
    % fun_plot_abs_yk(Yk_sample, using_simulation, dim)







    %% smoother plot
    
    vis_GA_x(bayes_first, xk_sample, dim, iter_first, strct_min_max, using_simulation, t_min);
    vis_GA_x(bayes_last, xk_sample, dim, iter_last, strct_min_max, using_simulation, t_min);




    %% GC
%     bayes_last = Bayes_Step{5};
    if using_simulation==0
        vis_GA_gc(bayes_first.GC_val, cfg_Init.GC_empirical, t_min, iter_first)
        vis_GA_gc(bayes_last.GC_val, cfg_Init.GC_empirical, t_min, iter_last)

%         bayes_1 = Bayes_Step{3};
%         vis_GA_x_gcohempr(bayes_1, xk_sample, dim, iter_last, strct_min_max, using_simulation, t_min, cfg_Init.GC_empirical);

    end

end



%% plot error tot iteration
% vis_GA_error_TOT(abs_error_TOT_iter, dim, Iter)
% close all

%% W plot
% vis_GA_W(W_sim, W_iter, dim, using_simulation, use_bias)
close all

%% mu plot
% vis_GA_mu(cnfg_sim.mu, MU_iter, dim, using_simulation)

%% cholseky plot

if using_simulation==1
    %% Calc Cholesky decomposition - plot
    vis_GA_cholesky(xk_sample, Yk, cnfg_sim, W_iter, MU_iter, dim, bayes_first, bayes_last, Iter)
    close all

elseif using_simulation==0

    w_type = 'init';
    x_tick_init = [-2, 0, 2];
    x_lim_init = [-4 4];
    y_tick_init = [-2, 0, 2];
    y_lim_init = [-4 4];
    W_init = squeeze(W_iter(1, :, :));
    x_inint = bayes_first.MU_smoother;
    mu_init = MU_iter(1, :).';
    L = cnfg_sim.L;
    vis_GA_cholesky_real(x_inint, Yk, L, W_init, mu_init, w_type, dim, x_tick_init, x_lim_init, y_tick_init, y_lim_init)
    close all


    w_type = 'final';
    x_tick_fin = [-2, 0, 2];
    x_lim_fin = [-4 4];
    y_tick_fin = [-2, 0, 2];
    y_lim_fin = [-4 4];
    x_final = bayes_last.MU_smoother;
    W_final = squeeze(W_iter(Iter, :, :));
    mu_final = MU_iter(Iter, :).';
    L = cnfg_sim.L;
    vis_GA_cholesky_real(x_final, Yk, L, W_final, mu_final, w_type, dim,  x_tick_fin, x_lim_fin, y_tick_fin, y_lim_fin)
    close all


end



%% nested functions
% load data
function [cfg_Init, Y_k]= gc_load_yk_cfg(ch_num, num_switch)

str_dict = '.\load_folder\';

switch ch_num
    case 32
        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr11_sr256_winSec8_segNum8_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr11_sr256_winSec8_segNum8_TRIAL', ch_num);

        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec32_segNum32_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec32_segNum32_TRIAL', ch_num);

    case 20
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr12_sr256_winSec40_segNum40_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum20_fr12_sr256_winSec40_segNum40_TRIAL', ch_num);

    case 10
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL_good', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL_good', ch_num);

        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);

    case 8
        str_name_1 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec16_segNum16_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec16_segNum16_TRIAL', ch_num);

    case 5
        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec8_segNum8_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec8_segNum8_TRIAL', ch_num);

        % 'Anethesia_switch_Yk_chNum5_fr12_sr250_winSec10_segNum10_TRIAL.mat';
        % 'Anethesia_switch_cfg_Init_chNum5_fr12_sr250_winSec10_segNum10_TRIAL.mat';

        str_name_1 = sprintf('Anethesia_switch_%d_Yk_chNum%d_fr12_sr250_winSec10_segNum10_TRIAL', num_switch, ch_num);
        str_name_2 = sprintf('Anethesia_switch_%d_cfg_Init_chNum%d_fr12_sr250_winSec10_segNum10_TRIAL', num_switch, ch_num);

        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr12_sr256_winSec10_segNum10_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr12_sr256_winSec10_segNum10_TRIAL', ch_num);

    case 3
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec15_segNum15_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec15_segNum15_TRIAL', ch_num);

        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec10_segNum10_TRIAL', ch_sim);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec10_segNum10_TRIAL', ch_sim);
end

str_load_1 = sprintf('%s%s', str_dict, str_name_1);
str_load_2 = sprintf('%s%s', str_dict, str_name_2);


load(str_load_1);
load(str_load_2);
% yk = load(str_load_1);
% Y_k = cell2mat(struct2cell(yk));
% b = load(str_load_2);
% b2 = struct2cell(b);
end

% generate initial values for W
function [W] = gc_gen_W_init(coef_w, ch_num)
w_1 = coef_w*(0 + 1*rand(ch_num, 1));
w_1(1,1) = 1;
w_2 = coef_w*(0 +1*randn(ch_num, 1));
w_2(1,1) = 1;
w_3 = coef_w*randn(ch_num, 1);
w_3(1,1) = 0;

% if use_bias ==0 && use_bias_2_end==0
%     w_3 = 0;
% end

W = zeros(ch_num, 2);
W(:, 1) = w_1;
W(:, 2) = w_2;
W(:, 3) = w_3;
save('W.txt', 'W', '-ascii')

end

% generate simulation data
function [xk_sample, Yk_sample, cnfg_sim] = gc_sim_data(dim, QQ_sim, L_init, KK, Y_k, ch_num, switch_num)
%% data simulation
cnfg_sim.dim = dim;
cnfg_sim.x0 = [0;0];
cnfg_sim.F = eye(2);
cnfg_sim.D = [0;0];
cnfg_sim.Q = QQ_sim;
cnfg_sim.L = L_init;
cnfg_sim.incongruent_vec = ones(KK, 1);
cnfg_sim.W = load('W.txt');
muu = mean(Y_k).';
cnfg_sim.mu = muu(1:ch_num, 1);

% call function to simulate data
[xk_sample, Yk_sample] = SIM_samp_x_yk_switch(cnfg_sim, switch_num);


str_save = sprintf('X_Y_sim%d.mat', dim);
save(str_save, 'xk_sample', 'Yk_sample', 'cnfg_sim')
end

% print info
function []= gc_print_fun_1(cfg_Init, ch_num)

win_sec_GCoh_empircal = cfg_Init.cfgGen.win_sec;
win_sec_observation = cfg_Init.cfgGen.win_sec / cfg_Init.cfgGen.seg_num;
fprintf('\n Number of channel: %d' , ch_num);
fprintf(' \n Length of Window for calc GCOh empirical: %d sec', win_sec_GCoh_empircal)
fprintf(' \n Length of subWindow for calc Cov mat of GCOh empirical (Window for SS-Coh): %d sec', win_sec_observation)
fprintf('\n------------------------------\n')

end

% generate initial vlaues parameters
function [init_Param, update_coef_W, W_sim] = gc_init_param(ParamaUpdate, cnfg_sim, coef_noise_mu, col_W_1, col_W_2, col_W_3, coef_noise_W, ch_num, dim)
count = 0;
%% Q
fprintf('\n Updating Parameters: \n')
init_Param = [];
if ParamaUpdate.sigma_update == 1
    count = count + 1;
    fprintf('%d- Sigma \n', count)

    init_Param.Q = 100*[1 0; 0 1].*10^-3;

elseif ParamaUpdate.sigma_update == 0
    init_Param.Q = cnfg_sim.Q;

end

if dim==1
    init_Param.Q = init_Param.Q(1, 1);
end

%% mu
if ParamaUpdate.mu_update == 1
    count = count + 1;
    fprintf('%d- MU \n', count)

    init_Param.mu = cnfg_sim.mu - coef_noise_mu.*(0.5 + 1i*0.5);

elseif ParamaUpdate.mu_update == 0
    init_Param.mu = cnfg_sim.mu;
end

%% L
if ParamaUpdate.L_update == 1
    count = count + 1;
    fprintf('%d- L \n', count)

    load('L1');
    init_Param.L  = L1;

elseif ParamaUpdate.L_update == 0
    init_Param.L  = cnfg_sim.L;
end

%% W
if ParamaUpdate.W_update == 1
    count = count + 1;
    fprintf('%d- W\n', count)

    W_sim = cnfg_sim.W;
    mw = size(W_sim);
    update_coef_W = zeros(mw);
    if exist('col_W_1', 'var')
        update_coef_W(col_W_1, 1) =1;
    end
    if exist('col_W_2', 'var')
        update_coef_W(col_W_2, 2) =1;
    end
    if exist('col_W_3', 'var')
        update_coef_W(col_W_3, 3) =1;
    end

    str_noise = sprintf('noise_%d_3.mat', ch_num);

    if isfile(str_noise)

        %         mu = 0; sigma = 1;
        %         pdN = makedist('Normal',mu,sigma);
        %         pdHN = truncate(pdN,0,inf);
        %         step = 1/ch_num;
        %         X = (4:step:5-step)';
        %         Xt = X(X>0);
        %         a_1 = pdf(pdHN,Xt);
        %         a_2 = pdf(pdHN,Xt);
        %         a_3 = pdf(pdHN,Xt);
        %         nn = [a_1, a_2, a_3];


        noise = randn(ch_num, 3);
        noise_source = noise;
        %         noise_source = load(str_noise);
    else
        noise = randn(ch_num, 3);
        str_save = sprintf('noise_%d_3.mat', ch_num);
        save(str_save, 'noise')
    end

    %     noise_source = cell2mat(struct2cell(noise_source));
    W_noise = coef_noise_W.*noise_source.*update_coef_W;

    W_inint = W_sim + W_noise;

    if dim==1
        W_inint(:, 2) = [];
        update_coef_W(:, 2) = [];
        W_sim(:, 2) = [];
    end

    init_Param.W = W_inint;

elseif ParamaUpdate.W_update == 0
    init_Param.W = cnfg_sim.W;
end
fprintf('------------------------------\n')

end

