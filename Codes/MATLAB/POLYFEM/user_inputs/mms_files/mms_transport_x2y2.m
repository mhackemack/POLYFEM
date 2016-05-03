function [code_rhs, code_sol, code_angsol] = mms_transport_x2y2
syms x y Lx Ly sigma_t mu eta
% Functional forms
ang_sol = x*y*(Lx-x)*(Ly-y);
s_sol = 2*pi*x*y*(Lx-x)*(Ly-y);
f_func = mu*diff(ang_sol,x) + eta*diff(ang_sol,y) + sigma_t*ang_sol;
% Generate Code Output
code_rhs = matlabFunction(f_func);
code_sol = matlabFunction(s_sol);
code_angsol = matlabFunction(ang_sol);