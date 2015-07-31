lin_func = @(x,y,z) 1-3*x+2*y+4*z;
% Make Cube Cell Info
v=[0,0,0;1,0,0;1,1,0;0,1,0;0,0,1;1,0,1;1,1,1;0,1,1];
f=cell(6,1);
f{1}=[1,2,3,4];
f{2}=[5,8,7,6];
f{3}=[1,5,6,2];
f{4}=[2,6,7,3];
f{5}=[3,7,8,4];
f{6}=[4,8,5,1];
% Make Tet Cell Info
vv = [0,0,0;1,0,0;0,1,0;0,0,1];
ff = cell(4,1);
ff{1} = [1,2,3];
ff{2} = [1,2,4];
ff{3} = [1,4,3];
ff{4} = [2,3,4];
% Get Cube Interp Info
[qx_c,qw_c,m_c]=PWLD_quad_gen(v,2,f);
vc_vals = lin_func(v(:,1),v(:,2),v(:,3));
qc_vals = lin_func(qx_c(:,1),qx_c(:,2),qx_c(:,3));
mc_vals = m_c*vc_vals;
M_c = PWLD_volume(v,f,[1,0,0]);
for i=1:length(f)
    [MM{i},GG{i}] = PWLD_surface_ind2(v,f{i});
end
% Get Tet Interp Info
[qx_t,qw_t,m_t]=PWLD_quad_gen(vv,2,ff);
vt_vals = lin_func(vv(:,1),vv(:,2),vv(:,3));
qt_vals = lin_func(qx_t(:,1),qx_t(:,2),qx_t(:,3));
mt_vals = m_t*vt_vals;
M_t = PWLD_volume(vv,ff,[1,0,0]);

% 2D Stuff
lin_func_2 = @(x,y) 1-3*x+2*y;
v2 = [0,0;1,0;1,1;0,1];
[qx,qw,m2] = PWLD_quad_gen(v2,2);
v2_vals = lin_func_2(v2(:,1),v2(:,2));
q2_vals = lin_func_2(qx(:,1),qx(:,2));
m2_vals = m2*v2_vals;
[M2,~,G2] = PWLD_volume(v2);
return