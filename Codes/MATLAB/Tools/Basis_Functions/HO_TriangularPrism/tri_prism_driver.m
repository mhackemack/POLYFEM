%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Higher-Order Triangular Prism Test Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tri_prism_driver()
clear; close all; clc;
% Reference user data
% ------------------------------------------------------------------------------
npts = 10;
eval_nodes = get_ref_eval_nodes(npts);
qx = eval_nodes(:,1);
qy = eval_nodes(:,2);
qz = eval_nodes(:,3);
% Form reference triangular prism data
% ------------------------------------------------------------------------------
tri_ref_verts = [0,0;1,0;0,1];
prism_ref_verts = [tri_ref_verts,zeros(3,1);tri_ref_verts,ones(3,1)];

% Check linear space
% ------------------------------------------------------------------------------
prism_nodes = get_prism_ref_nodes(1);
pnx = prism_nodes(:,1);
pny = prism_nodes(:,2);
pnz = prism_nodes(:,3);
lin_vals = eval_ref_prism_vals(1, eval_nodes);
fprintf('-> Verifying linear basis functions.\n')
fprintf('      c-constraint = %s!\n',pass_fail_string(1-sum(lin_vals,2)))
fprintf('      x-constraint = %s!\n',pass_fail_string(qx-lin_vals*pnx))
fprintf('      y-constraint = %s!\n',pass_fail_string(qy-lin_vals*pny))
fprintf('      z-constraint = %s!\n',pass_fail_string(qz-lin_vals*pnz))
fprintf('     xz-constraint = %s!\n',pass_fail_string(qx.*qz-lin_vals*(pnx.*pnz)))
fprintf('     yz-constraint = %s!\n',pass_fail_string(qy.*qz-lin_vals*(pny.*pnz)))
fprintf('     xy-constraint = %s!\n',pass_fail_string(qx.*qy-lin_vals*(pnx.*pny))) % This should fail!!!
% Check quadratic space
% ------------------------------------------------------------------------------
prism_nodes = get_prism_ref_nodes(2);
pnx = prism_nodes(:,1);
pny = prism_nodes(:,2);
pnz = prism_nodes(:,3);
quad_vals = eval_ref_prism_vals(2, eval_nodes);
fprintf('\n-> Verifying quadratic basis functions.\n')
fprintf('      c-constraint = %s!\n',pass_fail_string(1-sum(quad_vals,2)))
fprintf('      x-constraint = %s!\n',pass_fail_string(qx-quad_vals*pnx))
fprintf('      y-constraint = %s!\n',pass_fail_string(qy-quad_vals*pny))
fprintf('      z-constraint = %s!\n',pass_fail_string(qz-quad_vals*pnz))
fprintf('     xy-constraint = %s!\n',pass_fail_string(qx.*qy-quad_vals*(pnx.*pny)))
fprintf('     xz-constraint = %s!\n',pass_fail_string(qx.*qz-quad_vals*(pnx.*pnz)))
fprintf('     yz-constraint = %s!\n',pass_fail_string(qy.*qz-quad_vals*(pny.*pnz)))
fprintf('     xx-constraint = %s!\n',pass_fail_string(qx.*qx-quad_vals*(pnx.*pnx)))
fprintf('     yy-constraint = %s!\n',pass_fail_string(qy.*qy-quad_vals*(pny.*pny)))
fprintf('     zz-constraint = %s!\n',pass_fail_string(qz.*qz-quad_vals*(pnz.*pnz)))
fprintf('    xxz-constraint = %s!\n',pass_fail_string(qx.*qx.*qz-quad_vals*(pnx.*pnx.*pnz)))
fprintf('    yyz-constraint = %s!\n',pass_fail_string(qy.*qy.*qz-quad_vals*(pny.*pny.*pnz)))
fprintf('    xyz-constraint = %s!\n',pass_fail_string(qx.*qy.*qz-quad_vals*(pnx.*pny.*pnz)))
fprintf('    xzz-constraint = %s!\n',pass_fail_string(qx.*qz.*qz-quad_vals*(pnx.*pnz.*pnz)))
fprintf('    yzz-constraint = %s!\n',pass_fail_string(qy.*qz.*qz-quad_vals*(pny.*pnz.*pnz)))
fprintf('   xxzz-constraint = %s!\n',pass_fail_string(qx.*qx.*qz.*qz-quad_vals*(pnx.*pnx.*pnz.*pnz)))
fprintf('   yyzz-constraint = %s!\n',pass_fail_string(qy.*qy.*qz.*qz-quad_vals*(pny.*pny.*pnz.*pnz)))
fprintf('   xyzz-constraint = %s!\n',pass_fail_string(qx.*qy.*qz.*qz-quad_vals*(pnx.*pny.*pnz.*pnz)))
% Check cubic space
% ------------------------------------------------------------------------------
prism_nodes = get_prism_ref_nodes(3);
pnx = prism_nodes(:,1);
pny = prism_nodes(:,2);
pnz = prism_nodes(:,3);
cubic_vals = eval_ref_prism_vals(3, eval_nodes);
fprintf('\n-> Verifying cubic basis functions.\n')
fprintf('      c-constraint = %s!\n',pass_fail_string(1-sum(cubic_vals,2)))
fprintf('      x-constraint = %s!\n',pass_fail_string(qx-cubic_vals*pnx))
fprintf('      y-constraint = %s!\n',pass_fail_string(qy-cubic_vals*pny))
fprintf('      z-constraint = %s!\n',pass_fail_string(qz-cubic_vals*pnz))
fprintf('     xz-constraint = %s!\n',pass_fail_string(qx.*qz-cubic_vals*(pnx.*pnz)))
fprintf('     yz-constraint = %s!\n',pass_fail_string(qy.*qz-cubic_vals*(pny.*pnz)))
fprintf('     xx-constraint = %s!\n',pass_fail_string(qx.*qx-cubic_vals*(pnx.*pnx)))
fprintf('     yy-constraint = %s!\n',pass_fail_string(qy.*qy-cubic_vals*(pny.*pny)))
fprintf('     zz-constraint = %s!\n',pass_fail_string(qz.*qz-cubic_vals*(pnz.*pnz)))
fprintf('     xy-constraint = %s!\n',pass_fail_string(qx.*qy-cubic_vals*(pnx.*pny)))
fprintf('     xz-constraint = %s!\n',pass_fail_string(qx.*qz-cubic_vals*(pnx.*pnz)))
fprintf('     yz-constraint = %s!\n',pass_fail_string(qy.*qz-cubic_vals*(pny.*pnz)))
fprintf('     xx-constraint = %s!\n',pass_fail_string(qx.*qx-cubic_vals*(pnx.*pnx)))
fprintf('     yy-constraint = %s!\n',pass_fail_string(qy.*qy-cubic_vals*(pny.*pny)))
fprintf('     zz-constraint = %s!\n',pass_fail_string(qz.*qz-cubic_vals*(pnz.*pnz)))
fprintf('    xxz-constraint = %s!\n',pass_fail_string(qx.*qx.*qz-cubic_vals*(pnx.*pnx.*pnz)))
fprintf('    yyz-constraint = %s!\n',pass_fail_string(qy.*qy.*qz-cubic_vals*(pny.*pny.*pnz)))
fprintf('    xyz-constraint = %s!\n',pass_fail_string(qx.*qy.*qz-cubic_vals*(pnx.*pny.*pnz)))
fprintf('    xzz-constraint = %s!\n',pass_fail_string(qx.*qz.*qz-cubic_vals*(pnx.*pnz.*pnz)))
fprintf('    yzz-constraint = %s!\n',pass_fail_string(qy.*qz.*qz-cubic_vals*(pny.*pnz.*pnz)))
fprintf('   xxzz-constraint = %s!\n',pass_fail_string(qx.*qx.*qz.*qz-cubic_vals*(pnx.*pnx.*pnz.*pnz)))
fprintf('   yyzz-constraint = %s!\n',pass_fail_string(qy.*qy.*qz.*qz-cubic_vals*(pny.*pny.*pnz.*pnz)))
fprintf('   xyzz-constraint = %s!\n',pass_fail_string(qx.*qy.*qz.*qz-cubic_vals*(pnx.*pny.*pnz.*pnz)))
fprintf('    xxx-constraint = %s!\n',pass_fail_string(qx.*qx.*qx-cubic_vals*(pnx.*pnx.*pnx)))
fprintf('    yyy-constraint = %s!\n',pass_fail_string(qy.*qy.*qy-cubic_vals*(pny.*pny.*pny)))
fprintf('    zzz-constraint = %s!\n',pass_fail_string(qz.*qz.*qz-cubic_vals*(pnz.*pnz.*pnz)))
fprintf('    xxy-constraint = %s!\n',pass_fail_string(qx.*qx.*qy-cubic_vals*(pnx.*pnx.*pny)))
fprintf('    xyy-constraint = %s!\n',pass_fail_string(qx.*qy.*qy-cubic_vals*(pnx.*pny.*pny)))
% Check quartic space
% ------------------------------------------------------------------------------
prism_nodes = get_prism_ref_nodes(4);
pnx = prism_nodes(:,1);
pny = prism_nodes(:,2);
pnz = prism_nodes(:,3);
quartic_vals = eval_ref_prism_vals(4, eval_nodes);
fprintf('\n-> Verifying quartic basis functions.\n')
fprintf('      c-constraint = %s!\n',pass_fail_string(1-sum(quartic_vals,2)))
fprintf('      x-constraint = %s!\n',pass_fail_string(qx-quartic_vals*pnx))
fprintf('      y-constraint = %s!\n',pass_fail_string(qy-quartic_vals*pny))
fprintf('      z-constraint = %s!\n',pass_fail_string(qz-quartic_vals*pnz))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_eval_nodes(n)
out = zeros(n,3);
for i=1:n
    finished = false;
    while ~finished
        t = rand(1,3);
        if t(1) < 1 - t(2)
            out(i,:) = t;
            finished = true;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_prism_ref_nodes(deg)
tri_nodes = get_tri_ref_nodes(deg);
n = size(tri_nodes,1);
z = linspace(0,1,deg+1);
out = [];
for i=1:length(z)
    t = [tri_nodes,z(i)*ones(n,1)];
    out = [out;t];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_tri_ref_nodes(deg)
if deg == 1
    out = [0,0;1,0;0,1];
elseif deg == 2
    out = [0,0;1,0;0,1;.5,0;.5,.5;0,.5];
elseif deg == 3
    out = [0,0;1,0;0,1;1/3,0;2/3,0;2/3,1/3;1/3,2/3;0,2/3;0,1/3;1/3,1/3];
elseif deg == 4
    out = [0,0;1,0;0,1;1/4,0;2/4,0;3/4,0;3/4,1/4;2/4,2/4;1/4,3/4;0,3/4;0,2/4;0,1/4;1/4,1/4;2/4,1/4;1/4,2/4];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = eval_ref_prism_vals(deg, x)
tri_vals = tri2D_ref_vals(deg, x(:,1:2)); ntri = size(tri_vals,2);
z_vals = get_1D_values(deg, linspace(0,1,deg+1), x(:,3));
out = [];
for i=1:(deg+1)
    out = [out,tri_vals.*(z_vals(:,i)*ones(1,ntri))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = eval_ref_prism_grad(deg, x)
if deg == 1
    
elseif deg == 2
    
elseif deg == 3
    
elseif deg == 4
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = compute_tri_prism_jacobian(verts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = tri2D_ref_vals(deg, x)
s = x(:,1); t = x(:,2);
if deg == 1
    out = [1-s-t,s,t];
elseif deg == 2
    out = [  2*s.^2+4*s.*t+2*t.^2-3*s-3*t+1,...
             2*s.^2-s,...
             2*t.^2-t,...
             -4*s.^2-4*s.*t+4*s,...
             4*s.*t,...
             -4*t.^2-4*s.*t+4*t...
          ];
elseif deg == 3
    out = [9/2*(1-s-t).*(1/3-s-t).*(2/3-s-t),...
           9/2*s.*(1/3-s).*(2/3-s),...
           9/2*t.*(1/3-t).*(2/3-t),...
           27/2*s.*(1-s-t).*(2/3-s-t),...
          -27/2*s.*(1/3-s).*(1-s-t),...
          -27/2*s.*t.*(1/3-s),...
          -27/2*s.*t.*(1/3-t),...
          -27/2*t.*(1-s-t).*(1/3-t),...
           27/2*t.*(1-s-t).*(2/3-s-t),...
           27*s.*t.*(1-s-t)];
%     out = [(1-s-t).*(1-3*s-3*t).*(1-3/2*s-3/2*t),...
%             s.*(1-3*s).*(1-3/2*s),...
%             t.*(1-3*t).*(1-3/2*t),...
%             9*s.*(1-s-t).*(1-3/2*s-3/2*t),...
%             -9/2*s.*(1-s-t).*(1-3*s),...
%             -9/2*s.*t.*(1-3*s),...
%             -9/2*s.*t.*(1-3*t),...
%              9/2*(1-s-t).*t.*(1-3*s),...
%             9*(1-s-t).*t.*(1-3/2*s-3/2*t),...
%             27*s.*t.*(1-s-t)];
elseif deg == 4
    out = [64/6*(1-s-t).*(1/4-s-t).*(2/4-s-t).*(3/4-s-t),...
          -64/6*s.*(1/4-s).*(2/4-s).*(3/4-s),...
          -64/6*t.*(1/4-t).*(2/4-t).*(3/4-t),...
           256/6*s.*(1-s-t).*(2/4-s-t).*(3/4-s-t),...
          -768/12*s.*(1-s-t).*(1/4-s).*(3/4-s-t),...
           256/6*s.*(1-s-t).*(1/4-s).*(2/4-s),...
           256/6*s.*t.*(1/4-s).*(2/4-s),...
           768/12*s.*t.*(1/4-s).*(1/4-t),...
           256/6*s.*t.*(1/4-t).*(2/4-t),...
           256/6*t.*(1-s-t).*(1/4-t).*(2/4-t),...
          -768/12*t.*(1-s-t).*(1/4-t).*(3/4-s-t),...
           256/6*t.*(1-s-t).*(2/4-s-t).*(3/4-s-t),...
           768/6*s.*t.*(1-s-t).*(3/4-s-t),...
          -768/6*s.*t.*(1-s-t).*(1/4-s),...
          -768/6*s.*t.*(1-s-t).*(1/4-t)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = tri2D_ref_grad(deg, xx)
if deg == 1
    out = [-1,-1;1,0;0,1];
elseif deg == 2
    s = xx(:,1); t = xx(:,2);
    out = [  4*s+4*t-3,  4*s+4*t-3;...
             4*s-1,         0;...
             0,                4*t-1;...
             -8*s-4*t+4,-4*s;...
             4*t,           4*s;...
             -4*t,    -8*t-4*s+4];
elseif deg == 3
    s = xx(:,1); t = xx(:,2);
    out = [ - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 1/3) - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 2/3) - (9*(s + t - 1/3)*(s + t - 2/3))/2, - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 1/3) - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 2/3) - (9*(s + t - 1/3)*(s + t - 2/3))/2;...
                                                                  (9*(s - 1/3)*(s - 2/3))/2 + (9*s*(s - 1/3))/2 + (9*s*(s - 2/3))/2, 0;...
                                                                                                                                  0, (9*(t - 1/3)*(t - 2/3))/2 + (9*t*(t - 1/3))/2 + (9*t*(t - 2/3))/2;...
                                                   (27*s*(s + t - 1))/2 + (27*s*(s + t - 2/3))/2 + (27*(s + t - 1)*(s + t - 2/3))/2, (27*s*(s + t - 1))/2 + (27*s*(s + t - 2/3))/2;...
                                                         - (27*s*(s + t - 1))/2 - (27*(s - 1/3)*(s + t - 1))/2 - (27*s*(s - 1/3))/2,-(27*s*(s - 1/3))/2;...
                                                                                                    (27*s*t)/2 + (27*t*(s - 1/3))/2, (27*s*(s - 1/3))/2;...
                                                                                                                 (27*t*(t - 1/3))/2, (27*s*t)/2 + (27*s*(t - 1/3))/2;...
                                                                                                                -(27*t*(t - 1/3))/2, - (27*t*(s + t - 1))/2 - (27*(t - 1/3)*(s + t - 1))/2 - (27*t*(t - 1/3))/2;...
                                                                                      (27*t*(s + t - 1))/2 + (27*t*(s + t - 2/3))/2, (27*t*(s + t - 1))/2 + (27*t*(s + t - 2/3))/2 + (27*(s + t - 1)*(s + t - 2/3))/2;...
                                                                                                        - 27*t*(s + t - 1) - 27*s*t, - 27*s*(s + t - 1) - 27*s*t];
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_1D_values(ord, v, x)
if ord == 0
    out = ones(length(x),1);
elseif ord == 1
    out = [(v(2)-x), (x-v(1))]/(v(2) - v(1));
elseif ord == 2
    out = [(x-v(2)).*(x-v(3))/(v(1)-v(2))/(v(1)-v(3)),...
           (x-v(1)).*(x-v(3))/(v(2)-v(1))/(v(2)-v(3)),...
           (x-v(1)).*(x-v(2))/(v(3)-v(1))/(v(3)-v(2))];
elseif ord == 3
    out = [(x-v(2)).*(x-v(3)).*(x-v(4))/(v(1)-v(2))/(v(1)-v(3))/(v(1)-v(4)),...
           (x-v(1)).*(x-v(3)).*(x-v(4))/(v(2)-v(1))/(v(2)-v(3))/(v(2)-v(4)),...
           (x-v(1)).*(x-v(2)).*(x-v(4))/(v(3)-v(1))/(v(3)-v(2))/(v(3)-v(4)),...
           (x-v(1)).*(x-v(2)).*(x-v(3))/(v(4)-v(1))/(v(4)-v(2))/(v(4)-v(3))];
else
    out = ones(length(x),ord+1);
    for i=1:ord+1
        for j=1:ord+1
            if i==j, continue; end
            out(:,i) = out(:,i).*(x-v(j))/(v(i)-v(j));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = pass_fail_string(x)
if norm(x) < 1e-12
    out = 'PASSED';
else
    out = 'FAILED';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%