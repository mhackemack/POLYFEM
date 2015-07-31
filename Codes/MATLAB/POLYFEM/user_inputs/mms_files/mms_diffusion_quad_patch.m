function [code_rhs, code_sol] = mms_diffusion_quad_patch
syms x y z Lx Ly Lz c_diff sigma_a
% Functional forms
d_sol = 1-2*x+6*y+x^2-3*x*y+4*y^2;
f_func  = -1.0*c_diff*( diff(diff(d_sol,x),x) +  diff(diff(d_sol,y),y) );
% Generate Code Output
code_rhs = matlabFunction(f_func);
code_sol = matlabFunction(d_sol);