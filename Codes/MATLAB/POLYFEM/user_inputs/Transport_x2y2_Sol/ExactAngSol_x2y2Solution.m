function out = ExactAngSol_x2y2Solution(xx,~)
Lx = 1.0; Ly = 1.0;
x = xx(:,1); y = xx(:,2);
out = x.*y.*(Lx-x).*(Ly-y);