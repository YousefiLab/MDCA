function data_ana(T,ch1)
tt = T./(60*60);

figure
plot(tt , ch1);
xlabel('Time (hrs)')
xlim([0 tt(end)])
title('channel 1')

win_1 = ch1(1:256*8);
figure
plot((1:256*8) , win_1, 'LineWidth' , 2 , 'color' , [91, 13, 155]./255);
xlim([1 256*8])
xlabel('Samples')
title('Window 1')

seg_1 = win_1(1:128);
figure
plot((1:128),seg_1, 'LineWidth' , 3 , 'color' , [17, 104, 1]./255);
xlabel('Samples'), xlim([1 128])
title('Segment 1')


seg_2 = win_1(129:256);
figure
plot((129:256),seg_2, 'LineWidth' , 3 , 'color' , [17, 104, 1]./255);
xlabel('Samples'), xlim([129 256])
title('Segment 2')


seg_16 = win_1(2048-128 : 2048 );
figure
plot((2048-128 : 2048) , seg_16 , 'LineWidth' , 3 , 'color' , [17, 104, 1]./255);
xlim([2048-128  2048]);
xlabel('Samples')
title('Segment 16')

end