function [code_rhs, code_sol, code_ang_sol] = mms_transport_isotropic_quad_patch
syms x y Omegax Omegay t
% Functional forms
ang_sol = 1-2*x+6*y+x^2-3*x*y+4*y^2;
s_sol = int(ang_sol,t,0,2*pi);
f_func = Omegax*diff(ang_sol,x) + Omegay*diff(ang_sol,y);
% Generate Code Output
code_rhs = matlabFunction(f_func);
code_sol = matlabFunction(s_sol);
code_ang_sol = matlabFunction(ang_sol);