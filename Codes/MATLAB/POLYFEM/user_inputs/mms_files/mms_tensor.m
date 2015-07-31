function mms_tensor
clear all; clc; close all
% % syms x y Lx Ly x0 y0 c_diff sigma_a varia
% % 
% % U = 100*x*(Lx-x)*y*(Ly-y)*exp(-((x-x0)^2+(y-y0)^2)/varia)/(Lx*Ly)^2
% % 
% % forcing_fn = -1.0*c_diff*( diff(diff(U,x),x) +  diff(diff(U,y),y) ) + sigma_a *U
% % 
% % codeU = matlabFunction(U)
% % 
% % coderhs= matlabFunction(forcing_fn)

syms x y Lx Ly D1 D2 c_diff sigma_a K freq U grad_U flux forcing_fn2 forcing_fn

Lx_=1;
Ly_=Lx_;
freq_=3;
exact=@(x_,y_) sin(freq_*pi*x_/Lx_).*sin(freq*pi*y_/Ly_);

K=[ (x+1)^2+y^2 , -x*y ; -x*y, (x+1)^2]*c_diff
% K=eye(2)

U = sin(freq*pi*x/Lx).*sin(freq*pi*y/Ly);

grad_U = [diff(U,x); diff(U,y)]
flux = K * grad_U

forcing_fn =-1.0*( diff(flux(1),x) + diff(flux(2),y) ) + sigma_a *U

coderhs= matlabFunction(forcing_fn)

% forcing_fn2= -1.0*( diff(diff(U,x),x) +  diff(diff(U,y),y) ) + sigma_a *U
% forcing_fn2-forcing_fn

% xx=linspace(0,Lx_);
% yy=linspace(0,Ly_);
% [uu,vv]=meshgrid(xx,yy);
% zz=exact(uu,vv);
% figure(99)
% surf(uu,vv,zz)

return
end