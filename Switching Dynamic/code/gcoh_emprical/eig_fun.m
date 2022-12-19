function [eig_info , GC] = eig_fun(cross_spect)

m=size(cross_spect);
%%% m(1) = number of desired frequencies
%%% m(2) = number of windows

Eig_vec = cell(m(1),m(2));
Eig_val = cell(m(1),m(2));
GC = zeros(m(1),m(2));
for f_ind=1:m(1)
    
    for win_num=1:m(2)
    cross_spect_temp = (cross_spect{f_ind,win_num});
    [e_vector , e_value] = eig(cross_spect_temp);
    %%% e_vector = each column of this matrix is a eigenvector
    
    Eig_vec{f_ind,win_num} = e_vector;
    Eig_val{f_ind,win_num} = e_value;
    
    max_e_value = max(e_value(:));
    sum_e_value = sum(e_value(:));
    
    GC(f_ind,win_num) = max_e_value./sum_e_value;
    end
end
eig_info = cell(2,1);
eig_info{1,1} = Eig_val;
eig_info{2,1} = Eig_vec;

end