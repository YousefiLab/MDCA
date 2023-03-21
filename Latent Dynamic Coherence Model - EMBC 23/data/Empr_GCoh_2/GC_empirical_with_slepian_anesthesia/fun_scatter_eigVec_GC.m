function [] = fun_scatter_eigVec_GC(f_vec)
% this code calculate scatter plot for different angle of largest
% corresponding eigenvectors for step k and k+1 for different Global
% Cohernce value 
%
%Inputs
% f_vec = a vector that shows frequencies that we want to evalueate
%          (corresponding data of desired freq must be saved in current folder)
%
%
%
scatter_color = [174, 16, 232 ;
                244, 209, 66;
                46, 104, 50]./255;
            
scatter_lineWidth = [2 ;
                    2 ;
                    2];
            
f_num = length(f_vec);

ch_num = 64;
win_sec = 64;
seg_num = 64;
str_load_1 = sprintf('scatter_info_ch%d_win%d_seg%d_f' , ch_num, win_sec, seg_num);



figure
for f_ind = 1:f_num
    str_load_2 = sprintf('%d.mat' , f_vec(f_ind));
    str_load = strcat(str_load_1 , str_load_2);
    
    load(str_load)
    GC_vec(f_ind,:) = scatter_info.GC;
    angEigVec_vec(f_ind , :) = scatter_info.ang_eigVec;
    
    % plottig part
    sz = 40;
    s = scatter( angEigVec_vec(f_ind , :) , GC_vec(f_ind,:) , sz);
    s.LineWidth = scatter_lineWidth(f_ind , 1);
    s.MarkerEdgeColor = scatter_color(f_ind , :);
    s.MarkerFaceColor = scatter_color(f_ind , :);
    s.DisplayName = sprintf('f %d' , f_vec(f_ind));
    hold on
    
    
%     str_legen{f_ind} = sprintf('f %d' , f_vec(f_ind));
    
end
xlabel('Angle between Largest Corresponding Eigenvectors')
ylabel('Global Coherence')
xlim([0 pi./2])
ylim([0 1])

legend('show')


end