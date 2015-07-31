function mms_gaussian
clear all; clc; close all
syms x y z Lx Ly Lz x0 y0 z0 c_diff sigma_a varia

% U2 = 16*x*(Lx-x)*y*(Ly-y);
% U3 = 64*x*(Lx-x)*y*(Ly-y)*z*(Lz-z);
% UG2 = 100*x*(Lx-x)*y*(Ly-y)*exp(-((x-x0)^2+(y-y0)^2)/varia)/(Lx*Ly)^2;
% UG3 = 100*x*(Lx-x)*y*z*(Lz-z)*(Ly-y)*exp(-((x-x0)^2+(y-y0)^2+(z-z0)^2)/varia)/(Lx*Ly*Lz)^2;
% 
% forcing_fn_2  = -1.0*c_diff*( diff(diff(U2,x),x) +  diff(diff(U2,y),y) ) + sigma_a *U2;
% forcing_fn_3  = -1.0*c_diff*( diff(diff(U3,x),x) +  diff(diff(U3,y),y) ) + sigma_a *U3;
% forcing_fn_g2 = -1.0*c_diff*( diff(diff(UG2,x),x) +  diff(diff(UG2,y),y) ) + sigma_a *UG2;
% forcing_fn_g3 = -1.0*c_diff*( diff(diff(UG3,x),x) +  diff(diff(UG3,y),y) ) + sigma_a *UG3;
% 
% codeU2 = matlabFunction(U2)
% codeU3 = matlabFunction(U3)
% codeUG2 = matlabFunction(UG2)
% codeUG3 = matlabFunction(UG3)
% 
% coderhs2  = matlabFunction(forcing_fn_2)
% coderhs3  = matlabFunction(forcing_fn_3)
% coderhsg2 = matlabFunction(forcing_fn_g2)
% coderhsg3 = matlabFunction(forcing_fn_g3)


Lx = 1;
Ly = Lx;
Lz = Lx;
x0=3*Lx/4;y0=x0;z0=Lz/4;
varia=Lx^2/100;

exact2=@(x,y) 16*x.*(Lx-x).*y.*(Ly-y);
exact3=@(x,y,z) 16*x.*(Lx-x).*y.*(Ly-y).*z.*(Lz-z);
exactg2=@(x,y) 100*x.*(Lx-x).*y.*(Ly-y).*exp(-((x-x0).^2+(y-y0).^2)/varia)/(Lx*Ly)^2;
exactg3=@(x,y,z) 100*x.*(Lx-x).*y.*(Ly-y).*z.*(Lz-z).*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)/varia)/(Lx*Ly*Lz)^2;
xx=linspace(0,Lx);
yy=linspace(0,Ly);
zz=linspace(0,Lz);

[uu,vv]=meshgrid(xx,yy);
figure(98)
surf(uu,vv,exact2(uu,vv))
figure(99)
surf(uu,vv,exactg2(uu,vv))

[uu,vv,ww]=meshgrid(xx,yy,zz);
figure(100)
slice(uu,vv,ww,exact3(uu,vv,ww),.5,.5,[])
figure(101)
slice(uu,vv,ww,exactg3(uu,vv,ww),.75,.75,[])

return
end