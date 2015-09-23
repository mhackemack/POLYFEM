function [code_rhs, code_sol, code_angsol] = mms_transport_isotropic
syms x y sigma_t x0 y0 Lx Ly Omegax Omegay varia t
% Functional forms
ang_sol = 100/Lx^2/Ly^2*x*y*(Lx-x)*(Ly-y);
% ang_sol = 100/Lx^2/Ly^2*x*y*(Lx-x)*(Ly-y)*exp(-((x-x0)^2+(y-y0)^2)/varia);
s_sol = int(ang_sol,t,0,2*pi);
f_func = Omegax*diff(ang_sol,x) + Omegay*diff(ang_sol,y) + sigma_t*ang_sol;
% Generate Code Output
code_rhs = matlabFunction(f_func);
code_sol = matlabFunction(s_sol);
code_angsol = matlabFunction(ang_sol);