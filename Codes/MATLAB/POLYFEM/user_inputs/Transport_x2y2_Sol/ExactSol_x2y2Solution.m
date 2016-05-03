function out = ExactSol_x2y2Solution(xx,~)
Lx = 1.0; Ly = 1.0;
x = xx(:,1); y = xx(:,2);
out = x.*y.*pi.*(Lx-x).*(Ly-y).*2.0;