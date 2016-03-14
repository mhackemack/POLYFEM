function out = RHSFunc_QuadSolution(xx,dir)
x = xx(:,1); y = xx(:,2);
sigma_t = 1.0;
mu = dir(1); eta = dir(2);
a = 1;
b = 1;
c = 1;
d = 1;
e = 1;
f = 1;
out = eta.*(c+d.*x+f.*y.*2.0)+mu.*(b+d.*y+e.*x.*2.0)+sigma_t.*(a+b.*x+c.*y+e.*x.^2+f.*y.^2+d.*x.*y);
