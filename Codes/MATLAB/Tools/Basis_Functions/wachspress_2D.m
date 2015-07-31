function varargout = wachspress_2D(verts, xin, faces, varargin)
nout = nargout;
% check if any input
if nargin <= 1, error('No inputs specified.'); end
% get input information
[nv,dim] = size(verts);
nx = size(xin,1);
% get Wachpress interpolant values/grads
varargout{1} = get_values(verts, xin, nv, nx);
if nout == 2, varargout{2} = get_grads(verts, xin, varargout{1}, nv, nx, dim); end
return
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
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_single_value(verts, xin, nv)
vals = zeros(1,nv);
out = zeros(1,nv);
xx = verts - ones(nv,1)*xin;
xdim = size(xin,2);
for i=1:nv
    if norm(verts(i,:) - xin) < 1e-14
        out(i) = 1.0;
        return
    end
    if i==1
        v = [nv,i,i+1]';
    elseif i==nv
        v = [i-1,i,1]';
    else
        v = [i-1,i,i+1]';
    end
    if xdim == 2
        A1 = 0.5*det([ones(1,3);[xin;verts(v(1:2),:)]']);
        if abs(A1) < 1e-14 && dot(xx(v(1),:),xx(v(2),:)) < 0
            out(v(1)) = norm(xx(v(2),:));
            out(v(2)) = norm(xx(v(1),:));
            out = out / sum(out);
            return
        end
        A2 = 0.5*det([ones(1,3);[xin;verts(v(2:3),:)]']);
        if abs(A2) < 1e-14 && dot(xx(v(2),:),xx(v(3),:)) < 0
            out(v(2)) = norm(xx(v(3),:));
            out(v(3)) = norm(xx(v(2),:));
            out = out / sum(out);
            return
        end
        Av = 0.5*det([ones(1,3);verts(v,:)']);
    else
        xd1 = [xin;verts(v(1:2),:)];
        xd2 = [xin;verts(v(2:3),:)];
        A1 = polygonArea3d(xd1);
        if abs(A1) < 1e-14 && dot(xx(v(1),:),xx(v(2),:)) < 0
            out(v(1)) = norm(xx(v(2),:));
            out(v(2)) = norm(xx(v(1),:));
            out = out / sum(out);
            return
        end
        A2 = polygonArea3d(xd2);
        if abs(A2) < 1e-14 && dot(xx(v(2),:),xx(v(3),:)) < 0
            out(v(2)) = norm(xx(v(3),:));
            out(v(3)) = norm(xx(v(2),:));
            out = out / sum(out);
            return
        end
        Av = polygonArea3d(verts(v,:));
    end
    vals(i) = Av / (A1*A2);
end
out = vals / sum(vals);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grads = get_grads(verts, xin, vals, nv, nx, dim)
grads = zeros(nv, dim, nx);
n = get_edge_normals(verts,nv);
h = get_vert_distances(verts,xin,n,nv,nx);
for j=1:nx
    grads(:,:,j) = get_single_grad(verts, xin(j,:), vals(j,:), nv, n, h(j,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = get_single_grad(verts, xin, val, nv, n, h)
grad = zeros(nv,2);
xx = verts - ones(nv,1)*xin;
R = zeros(nv,2);
for i=1:nv
    if norm(xx(i,:)) < 1e-14
        grad = nan;
        return
    end
end
for i=1:nv
    if i==1
        v = [nv,i];
    else
        v = [i-1,i];
    end
    if abs(h(v(1))) < 1e-14
        dir = verts(v(2),:) - verts(v(1),:) / norm(verts(v(2),:) - verts(v(1),:));
        grad(v,:) = [-dir;dir];
        return
    end
    R(i,:) = n(v(1),:)/h(v(1)) + n(v(2),:)/h(v(2));
end
for i=1:nv
    grad(i,:) = val(i)*(R(i,:) - val*R);
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_edge_normals(verts,nv)
out = zeros(size(verts,1), size(verts,2));
vdim = size(verts,2);
vc = mean(verts);
if vdim == 2
    for i=1:nv
        if i==nv
            v = [i,1];
        else
            v = [i,i+1];
        end
        vv = verts(v,:);
        dd = vv(2,:) - vv(1,:);
        out(i,:) = [-dd(2),dd(1)];
        if dot(vc,out(i,:)) < 0
            out(i,:) = -1.0*out(i,:);
        end
    end
else
    out = polygon3dNormalAngle(verts,1:nv);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_vert_distances(verts,xin,n,nv,nx)
out = zeros(size(xin,1),size(verts,1));
for j=1:nx
    for i=1:nv
        out(j,i) = dot(n(i,:),verts(i,:) - xin(j,:));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
