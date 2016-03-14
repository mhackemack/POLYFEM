function out = ExactSol_QuadSolution(xx,~)
x = xx(:,1); y = xx(:,2);
a = 1;
b = 1;
c = 1;
d = 1;
e = 1;
f = 1;
out = pi.*(a+b.*x+c.*y+e.*x.^2+f.*y.^2+d.*x.*y).*2.0;