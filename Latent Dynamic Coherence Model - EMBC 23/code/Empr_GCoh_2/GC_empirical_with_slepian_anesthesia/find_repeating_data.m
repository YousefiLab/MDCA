clc
clearvars -except data_det T Fs mat_repeat
close all



%% data
%%% loading data
if exist('data_det') ~=1
    load('eeganes07laplac250_detrend_all.mat');
    
    %% finding repeating data
    mat_repeat = zeros(64,64);
    for i=1:64
        data_1=  data_det(:,i);
        i
        for j=i:64
            data_2=  data_det(:,j);
            mat_repeat(i, j) = sum(abs(data_1- data_2));
        end
    end
    
end

save('mat_abs_dif.mat', 'mat_repeat');
for i=1:64
    
    vec_i = mat_repeat(i, i+1:64);
    min_vec = min(vec_i)
    max_vec = max(vec_i)
    
end

ind_0 = find(mat_repeat);
