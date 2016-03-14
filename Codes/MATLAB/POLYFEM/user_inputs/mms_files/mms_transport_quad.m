function [code_rhs, code_sol, code_angsol] = mms_transport_quad()
syms a b c d e f mu eta x y sigma_t
% Functional forms
ang_sol = a + b*x + c*y + d*x*y + e*x*x + f*y*y;
s_sol = 2*pi*(a + b*x + c*y + d*x*y + e*x*x + f*y*y);
f_func = mu*diff(ang_sol,x) + eta*diff(ang_sol,y) + sigma_t*ang_sol;
% Generate Code Output
code_rhs = matlabFunction(f_func);
code_sol = matlabFunction(s_sol);
code_angsol = matlabFunction(ang_sol);