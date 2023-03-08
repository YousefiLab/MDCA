function [] = ang_diff_phase_ev(sorted_eig_info , config, fnt_size)

eVal = sorted_eig_info{1,1}; 
eVec = sorted_eig_info{2,1}; 

m = size(eVec);
mm = size(eVec{1,1});

Phase_diff_Large = zeros(mm(1), m(2));
Phase_diff_SecLarge = zeros(mm(1), m(2));


for i=1:length(eVec)
    temp_eVec = eVec{1, i};
    vec_Larg_eVec = temp_eVec(:, 1);
    vec_SecLarg_eVec = temp_eVec(:, 2);
    
    % between [-pi pi]
    theta_Larg_eVec = angle(vec_Larg_eVec);
    theta_SecLarg_eVec = angle(vec_SecLarg_eVec);
    
    % Reference CPz
    theta_Larg_eVec = theta_Larg_eVec-theta_Larg_eVec(7);
    theta_SecLarg_eVec = theta_SecLarg_eVec-theta_SecLarg_eVec(7);
    
    Phase_diff_Large(:, i) = theta_Larg_eVec;
    Phase_diff_SecLarge(:, i) = theta_SecLarg_eVec;
    
end

struct_ph = [];

struct_ph.reference = 'CPz';
struct_ph.num_ch = 32;
struct_ph.chs = [1 2 3 4 5 6 7 11 12 13 14 17 18 19 24 25 29 33 34 35 36 37 38 40 41 42 43 48 49 54 55 60];
struct_ph.Phase_diff_Large = Phase_diff_Large;
struct_ph.Phase_diff_SecLarge = Phase_diff_SecLarge;

save('strct_phase_diff_fr12_ch32.mat', 'struct_ph');

end