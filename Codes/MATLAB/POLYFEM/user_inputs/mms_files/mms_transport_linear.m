function [code_rhs, code_sol, code_angsol] = mms_transport_linear()
syms a b c d e f mu eta x y sigma_t
% Functional forms
ang_sol = a*x + b*y + c*mu + d*eta + e;
s_sol   = 2*pi*(a*x + b*y + e);
f_func = mu*diff(ang_sol,x) + eta*diff(ang_sol,y) + sigma_t*ang_sol;
% Generate Code Output
code_rhs = matlabFunction(f_func);
code_sol = matlabFunction(s_sol);
code_angsol = matlabFunction(ang_sol);