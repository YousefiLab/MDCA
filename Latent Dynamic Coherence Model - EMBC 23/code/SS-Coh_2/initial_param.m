function [cfg_Init] = initial_param(Yk, y_interval, ch_rand, ab_mood, b_coef)

% ab_mood = a -> a will be estimated and b =1
% ab_mood = b -> b will be estimated and a =0

Y_k = zeros(y_interval(2)-y_interval(1)+1, length(ch_rand));
for i=1: length(ch_rand)
    
    Y_k(:, i) = Yk(y_interval(1):y_interval(2), ch_rand(i));
end

m = size(Y_k);

mu_init = mean(Y_k).';

%%  cross spectral
cross_spect_mat = zeros(m(2), m(2));
K = m(1);
for i=1:m(2)
    for j=1:m(2)
        
        X_i = Y_k(:,i);
        X_j = Y_k(:,j);
        
        cross_spect_mat(i,j) = (1/K)*sum(X_i.*conj(X_j));
    end
end

%% L, D
[L , D] = eig(cross_spect_mat);

% sort
D_sort = zeros(m(2));
L_sort = zeros(m(2));
for i=1 : m(2)
    D_sort(i , i) = D(end - (i-1) , end - (i-1));
    L_sort(: , i) = L(: , end - (i-1));
end
D_init = D_sort;
L_init = L_sort;

%%  estimate a and b


switch ab_mood
    case 'a'
        x0 = log(D_init(1,1));

        b_init = b_coef*ones(m(2), 1);
        a_init = zeros(m(2), 1);
        
        for i=2:m(2)
            a_init(i, 1) = log(D_init(i,i))-x0;
        end
        
    case 'new'
        x0 = 0;

        a_init = zeros(m(2), 1);
        b_init = ones(m(2), 1);
        for i=1:m(2)
            a_init(i, 1) = log(D_init(i,i));
        end
%         a_init(1:2) = a_init(1:2)-2; 
        
        
    case 'b'
        x0 = log(D_init(1,1));

        b_init = zeros(m(2), 1);
        a_init = zeros(m(2), 1);
        
        b_init(1, 1) =  1;
        for i=2:m(2)
            b_init(i, 1) = log(D_init(i,i))./x0;
        end
        
    case 'c'
        x0 = log(D_init(1,1));

        b_init = zeros(m(2), 1);
        a_init = 0.5*ones(m(2), 1);
        
        a_init(1,1) = 0;
        b_init(1, 1) =  1;
        for i=2:m(2)
            b_init(i, 1) = (log(D_init(i,i))-a_init(i, 1))./x0;
        end
end




cfg_Init.mu = mu_init;
cfg_Init.L = L_init;
cfg_Init.D = D_init;
cfg_Init.a = a_init;
cfg_Init.b = b_init;
cfg_Init.Y_k = Y_k;
cfg_Init.x0 = x0;

end