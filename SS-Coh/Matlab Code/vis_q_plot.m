function [] = vis_q_plot(cnfg_sim, q, dim)

color_est_mat = [ 181, 7, 80; 13, 93, 158]./255;



if dim==1
    q_sim_1 = cnfg_sim.Q(1,1);
    
    % call plot function
    str_tit = sprintf('Q update || starting piont %.4f | last point %.4f | real vlaue %.4f', q(1, 1) ,q(1, end) , q_sim_1);
    fun_nstd_one_vec(q_sim_1, q, str_tit, color_est_mat(1, :));
    
elseif dim==2
    m = size(q);
    
    for i=1: m(1)
        q_temp = squeeze(q(i, :, :));
        
        q_1(1, i) = q_temp(1, 1);
        q_2(1, i) = q_temp(2, 2);
    end
    
    q_sim_1 = cnfg_sim.Q(1,1);
    q_sim_2 = cnfg_sim.Q(2,2);
    
    % call plot function
    str_tit = sprintf('Q(1, 1) update || starting piont %.4f | last point %.4f | real vlaue %.4f', q_1(1, 1) ,q_1(1, end) , q_sim_1);
    fun_nstd_one_vec(q_sim_1, q_1, str_tit, color_est_mat(1, :));
    
    str_tit = sprintf('Q(2, 2) update || starting piont %.4f | last point %.4f | real vlaue %.4f', q_2(1, 1) ,q_2(1, end) , q_sim_2);
    fun_nstd_one_vec(q_sim_2, q_2, str_tit, color_est_mat(2, :));
end

end

%% neseted
% plot one vector

function [] = fun_nstd_one_vec(q_sim, q_vec, str_tit, color_vec)
Iter = length(q_vec);

figure
plot(1: Iter, q_vec, '--', 'LineWidth', 2, 'Color', color_vec);
hold on
plot(1: Iter, q_sim.*ones(1,Iter), 'LineWidth', 2, 'Color', [0, 0, 0]./255);
hold off

max_1 = max([q_vec, q_sim]);
min_1 = min([q_vec, q_sim]);

ylim([min_1-.05 max_1+.05])
xlim([1 Iter])

xlabel('iter'),

title(str_tit)
legend('estimated values', 'original value')

end

