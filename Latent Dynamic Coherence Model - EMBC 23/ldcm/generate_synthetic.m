function out = generate_synthetic(T)

%%%% centered freq of different channels
%%%% each channel has 3 part with different frequencies
% % % % % % % f_1 = [14 22 30];
% % % % % % % f_2 = [30 21 60];
% % % % % % % f_3 = [31 22 10];
% % % % % % % f_4 = [7 23 90];
% % % % % % % f_5 = [14 22 70];
% % % % % % % f_6 = [16 21 17];
% % % % % % % f_7 = [29 22 3];
% % % % % % % f_8 = [70 18 50];
% % % % % % % 
% % % % % % % %%%% coef specifies amplitude of each channel (each channel has 3 parts and finally 3 amplitudes)
% % % % % % % coef_1 = [1 2 .2];
% % % % % % % coef_2 = [1 1 .3];
% % % % % % % coef_3 = [1 2 .5];
% % % % % % % coef_4 = [.8 2 .5];
% % % % % % % coef_5 = [.9 1 .4];
% % % % % % % coef_6 = [.8 2 .8];
% % % % % % % coef_7 = [1 1 .6];
% % % % % % % coef_8 = [1 1 .4];

f_1 = [14 20 30];
f_2 = [14 10 30];
f_3 = [31 22 10];
f_4 = [7 23 90];
f_5 = [14 22 70];
f_6 = [16 21 17];
f_7 = [29 22 3];
f_8 = [70 18 50];

%%%% coef specifies amplitude of each channel (each channel has 3 parts and finally 3 amplitudes)
coef_1 = [1 .2 1];
coef_2 = [1 .2 1];
coef_3 = [1 2 .5];
coef_4 = [1 2 .5];
coef_5 = [1 1 .4];
coef_6 = [1 2 .8];
coef_7 = [1 1 .6];
coef_8 = [1 1 .4];

%%% noise_coef specifies the amplitude of additive Gaussian noise
noise_coef_1 = .5;
noise_coef_2 = .5;
noise_coef_3 = 4.0;
noise_coef_4 = 4.0;
noise_coef_5 = 4.0;
noise_coef_6 = 4.0;
noise_coef_7 = 4.0;
noise_coef_8 = 4.0;

t = T(1 , 1:floor(length(T)/3));

%%% Generating channels 
x_1 = [coef_1(1)*sin(2*pi*f_1(1)*t) coef_1(2)*sin(2*pi*f_1(2)*t)  coef_1(3)*sin(2*pi*f_1(3)*t) 0 0 ];
x_2 = [coef_2(1)*sin(2*pi*f_2(1)*t) coef_2(2)*sin(2*pi*f_2(2)*t)  coef_2(3)*sin(2*pi*f_2(3)*t) 0 0];
x_3 = [coef_3(1)*sin(2*pi*f_3(1)*t) coef_3(2)*sin(2*pi*f_3(2)*t)  coef_3(3)*sin(2*pi*f_3(3)*t) 0 0 ];
x_4 = [coef_4(1)*sin(2*pi*f_4(1)*t) coef_4(2)*sin(2*pi*f_4(2)*t)  coef_4(3)*sin(2*pi*f_4(3)*t) 0 0];
x_5 = [coef_5(1)*sin(2*pi*f_5(1)*t) coef_5(2)*sin(2*pi*f_5(2)*t)  coef_5(3)*sin(2*pi*f_5(3)*t) 0 0 ];
x_6 = [coef_6(1)*sin(2*pi*f_6(1)*t) coef_6(2)*sin(2*pi*f_6(2)*t)  coef_6(3)*sin(2*pi*f_6(3)*t) 0 0];
x_7 = [coef_7(1)*sin(2*pi*f_7(1)*t) coef_7(2)*sin(2*pi*f_7(2)*t)  coef_7(3)*sin(2*pi*f_7(3)*t) 0 0 ];
x_8 = [coef_8(1)*sin(2*pi*f_8(1)*t) coef_8(2)*sin(2*pi*f_8(2)*t)  coef_8(3)*sin(2*pi*f_8(3)*t) 0 0];

%%%% Adding noise to signals
data_syn(:,1) = x_1 + noise_coef_1*randn(1,length(x_1));
data_syn(:,2) = x_2 + noise_coef_2*randn(1,length(x_2));
data_syn(:,3) = x_3 + noise_coef_3*randn(1,length(x_3));
data_syn(:,4) = x_4 + noise_coef_4*randn(1,length(x_4));
data_syn(:,5) = x_5 + noise_coef_5*randn(1,length(x_5));
data_syn(:,6) = x_6 + noise_coef_6*randn(1,length(x_6));
data_syn(:,7) = x_7 + noise_coef_7*randn(1,length(x_7));
data_syn(:,8) = x_8 + noise_coef_8*randn(1,length(x_8));

%%% time plotting
% % % % % % % % % % % figure,
% % % % % % % % % % % subplot(421),plot(T,data_syn(:,1));
% % % % % % % % % % % % xlabel('time(sec)'),
% % % % % % % % % % % ylabel('Amp'),
% % % % % % % % % % % ylim([-25 25]),title('Channel 1');
% % % % % % % % % % % 
% % % % % % % % % % % subplot(422),plot(T,data_syn(:,2));
% % % % % % % % % % % % xlabel('time(sec)'),
% % % % % % % % % % % ylabel('Amp'),
% % % % % % % % % % % ylim([-25 25]),title('Channel 2');
% % % % % % % % % % % 
% % % % % % % % % % % subplot(423),plot(T,data_syn(:,3));
% % % % % % % % % % % % xlabel('time(sec)'),
% % % % % % % % % % % ylabel('Amp'),
% % % % % % % % % % % ylim([-25 25]),title('Channel 3');
% % % % % % % % % % % 
% % % % % % % % % % % subplot(424),plot(T,data_syn(:,4));
% % % % % % % % % % % % xlabel('time(sec)'),
% % % % % % % % % % % ylabel('Amp'),
% % % % % % % % % % % ylim([-25 25]),title('Channel 4');
% % % % % % % % % % % 
% % % % % % % % % % % subplot(425),plot(T,data_syn(:,5));
% % % % % % % % % % % % xlabel('time(sec)'),
% % % % % % % % % % % ylabel('Amp'),
% % % % % % % % % % % ylim([-25 25]),title('Channel 5');
% % % % % % % % % % % 
% % % % % % % % % % % subplot(426),plot(T,data_syn(:,6));
% % % % % % % % % % % % xlabel('time(sec)'),
% % % % % % % % % % % ylabel('Amp'),
% % % % % % % % % % % ylim([-25 25]),title('Channel 6');
% % % % % % % % % % % 
% % % % % % % % % % % subplot(427),plot(T,data_syn(:,7));
% % % % % % % % % % % xlabel('time(sec)'),ylabel('Amp'),ylim([-25 25]),title('Channel 7');
% % % % % % % % % % % 
% % % % % % % % % % % subplot(428),plot(T,data_syn(:,8));
% % % % % % % % % % % xlabel('time(sec)'),ylabel('Amp'),ylim([-25 25]),title('Channel 8');
% % % % % % % % % % % 
% % % % % % % % % % % suptitle('time plotting of 8 channels')
% % % % % % % % % % % % close

out = data_syn;
end