clc 
clear 
close all 

rng default  % For reproducibility
p1 = pearsrnd(0,1,-1,4,1000,1);
p2 = pearsrnd(0,1,0.75,3,1000,1);

figure
scatterhist(p1,p2)

k=1;