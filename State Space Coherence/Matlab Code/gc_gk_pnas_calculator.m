function [gc_pnas] = gc_gk_pnas_calculator(Y_k, win, over_lap)
% GC_GK_PNAS_CALCULATOR calculates Global Coherence measure based PNAS
% method for simulated data produced with the proposed model
% 
% INPUTS
% Y_k      : simulated data
% win      : number of steps that will be used to calculate GC
% over_lap : number of overlapped steps that will be used for calculate GC
% 
% OUTPUTS
% gc_pnas     :

m = size(Y_k);
% m(1) = number of points
% m(2) = number of channels (dimension)

if nargin==0
    error('the first input "Y_k" should be passed')
end

if nargin==1
    win = 2*m(2);
    over_lap = floor(.25*win);
end
if nargin==2
    over_lap = floor(.25*win);
end


% length of jumping with considering overlap
j_s = win - over_lap;  % Jumping Step

K = m(1); % number of points

% calculate g_k for each 
for i=1:floor((K-win)/j_s)
    
    
    for j=1:m(2) % number of channels
        vec_fft(j,:) = Y_k((i-1)*j_s + 1:(i-1)*j_s + win , j); % vector of fft at specifeid frequncy        
    end
    
    
    %% calculate cross spectral matrix 
    
    cross_spct = zeros(m(2), m(2));
    % a square matrix with length of number of channels
    
    for k=1:m(2)
        for l=1:m(2)
            cross_spct(k,l) = (vec_fft(k,:)*vec_fft(l,:)')/win;
        end
    end
    
    %% calculate eigeninfo (eigenvector [L], eigenvalue [D])
    [L , D] = eig(cross_spct); % decompose cross spectral matrix to eigenvalue and eigenvector matrix
    
    %% caclculate GC measure
    lambda = diag(D); % diagonal elements of eigenvalue matrix
    
    gc_pnas(i,1) = max(lambda)/sum(lambda); % Global Coherence measure based PNAS paper
end


end