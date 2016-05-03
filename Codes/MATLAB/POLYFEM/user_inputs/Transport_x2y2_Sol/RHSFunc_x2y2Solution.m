function out = RHSFunc_x2y2Solution(xx,dir)
Lx = 1.0; Ly = 1.0;
x = xx(:,1); y = xx(:,2);
sigma_t = 1.0;
mu = dir(1); eta = dir(2);
out = -eta.*(x.*y.*(Lx-x)-x.*(Lx-x).*(Ly-y))-mu.*(x.*y.*(Ly-y)-y.*(Lx-x).*(Ly-y))+sigma_t.*x.*y.*(Lx-x).*(Ly-y);
