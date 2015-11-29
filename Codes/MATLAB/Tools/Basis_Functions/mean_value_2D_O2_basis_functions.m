function varargout = mean_value_2D_O2_basis_functions(verts, xin, faces, varargin)
nout = nargout;
% check if any input
if nargin <= 1, error('No inputs specified.'); end
% get input information
[nv,dim] = size(verts);
nx = size(xin,1);
% get MV interpolant values/grads
varargout{1} = get_values(verts, xin, nv, nx);
if nout == 2
    varargout{2} = get_grads(verts, xin, varargout{1}, nv, nx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Function Lists
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vals = get_values(verts, xin, nv, nx)
vals = zeros(nx, nv);
for j=1:nx
    vals(j,:) = get_single_value(verts, xin(j,:), nv);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_single_value(verts, xin, nv)
out = zeros(1,nv);
xx = verts - ones(nv,1)*xin;
A = zeros(nv,1);
D = zeros(nv,1);
t = zeros(nv,1);
vdim = size(verts,2);
for i=1:nv
    if norm(xx(i,:)) < 1e-14
        out(i) = 1;
        return
    end
    if i==nv
        v = [i,1];
    else
        v = [i,i+1];
    end
    vt = [xx(v(1),:);xx(v(2),:)]';
    if vdim == 3
        vt = [vt,[0;0;0]];
    end
    A(i) = (det(vt));
    D(i) = dot(xx(v(1),:),xx(v(2),:));
    if abs(A(i)) < 1e-14 && D(i) < 0
        out(v(1)) = norm(xx(v(2),:));
        out(v(2)) = norm(xx(v(1),:));
        out = out / sum(out);
        return
    end
    t(i) = A(i) / (norm(xx(v(1),:))*norm(xx(v(2),:)) + D(i));
end
for i=1:nv
    if i==1
        v = [nv,1];
    else
        v = [i-1,i];
    end
    out(i) = sum(t(v)) / norm(xx(v(2),:));
end
out = out / sum(out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grads = get_grads(verts, xin, vals, nv, nx)
grads = zeros(nv, 2, nx);
for j=1:nx
    grads(:,:,j) = get_single_grad(verts, xin(j,:), vals(j,:), nv);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = get_single_grad(verts, xin, val, nv)
grad = zeros(nv,2);
xx = verts - ones(nv,1)*xin;
dim = length(xin);
R = zeros(nv,dim);
for i=1:nv
    if norm(xx(i,:)) < 1e-14
        grad = nan;
        return
    end
    if i==1
        v = [nv,i,i+1];
    elseif i==nv
        v = [i-1,i,1];
    else
        v = [i-1,i,i+1];
    end
    r1 = norm(xx(v(1),:));
    r2 = norm(xx(v(2),:));
    r3 = norm(xx(v(3),:));
    ang1 = acos(dot(xx(v(2),:),xx(v(1),:))/(r2*r1));
    
    if abs(sin(ang1)) < 1e-14
        dir = verts(v(2),:) - verts(v(1),:) / norm(verts(v(2),:) - verts(v(1),:));
        grad(v(1:2),:) = [-dir;dir];
        return
    end
    
	ang2 = acos(dot(xx(v(3),:),xx(v(2),:))/(r3*r2));
    e1 = xx(v(1),:) / r1;
    e2 = xx(v(2),:) / r2;
    e3 = xx(v(3),:) / r3;
    t1 = tan(ang1/2);
    t2 = tan(ang2/2);
    c1 = e1/r1 - e2/r2;
    c2 = e2/r2 - e3/r3;
    if dim == 2
        R(i,:) = t1/(t1+t2)*[-c1(2),c1(1)]/sin(ang1) + t2/(t1+t2)*[-c2(2),c2(1)]/sin(ang2) + e2/r2;
    elseif dim == 3
        
%         R(i,:) = (t1/(t1+t2))/sin(ang1)*[-c1(2),c1(1)] + (t2/(t1+t2))/sin(ang2)*[-c2(2),c2(1)] + e2/r2;
    end
end
for i=1:nv
    grad(i,:) = val(i)*(R(i,:) - val*R);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

