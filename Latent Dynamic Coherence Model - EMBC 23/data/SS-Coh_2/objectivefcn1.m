function f = objectivefcn1(W, MU_smoother, SIGMA_smoother, Yk, mu, L, incongruent_vec, dim)

f = 0;

for k = -10:10
    
    f = f + exp(-(x(1)-x(2))^2 - 2*x(1)^2)*cos(x(2))*sin(2*x(2));
end