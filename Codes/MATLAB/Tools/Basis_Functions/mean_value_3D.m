function varargout = mean_value_3D(verts, xin, faces, varargin)
% Get input/output information
% ----------------------------
if nargin <= 1, error('No inputs specified.'); end
nv = size(verts,1);
nx = size(xin,1);
dim1 = size(verts,2);
dim2 = size(xin,2);
if dim1 ~= dim2, error('Dimensions between inputs does not match.'); end

% Get MV interpolant values/grads
% -------------------------------
p = get_face_planes(verts, faces);
if nargout == 1
    varargout{1} = get_values(verts, faces, p, xin, nv, nx);
elseif nargout == 2
    [varargout{1},varargout{2}] = get_values(verts, faces, p, xin, nv, nx, dim1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Function Lists
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = get_values(verts, faces, planes, xin, nv, nx, dim)
nout = nargout;
vals = zeros(nx, nv);
if nout == 2
    grads = zeros(nx, nv, dim);
end
if nout == 1
    for j=1:nx
        vals(j,:) = get_single_value(verts, faces, planes, xin(j,:), nv);
    end
else
    for j=1:nx
        [vals(j,:), grads(j,:,:)] = get_single_value(verts, faces, planes, xin(j,:), nv);
    end
end
% outputs
varargout{1} = vals;
if nout == 2, varargout{2} = grads; end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = get_single_value(verts, faces, p, xin, nv)
nout = nargout;
if nout == 2
    R = zeros(nv,3);
    grads = zeros(nv,3);
end
w = zeros(1,nv);
xx = verts - ones(nv,1)*xin;
u = zeros(nv,3);
d = zeros(nv,1);
for i=1:nv
    u(i,:) = xx(i,:);
    d(i) = norm(u(i,:));
    if d(i) < 1e-14
        w(i) = 1;
        varargout{1} = w;
        if nout==2
            grads(:,:) = nan;
            varargout{2} = grads;
        end
        return
    end
    u(i,:) = u(i,:) / d(i);
end
for f=1:length(faces)
    face = faces{f};
    if abs(distancePointPlane(xin,p{f})) < 1e-13
        if nout == 1
            w(1,face) = mean_value_2D(verts(face,:), xin);
            varargout{1} = w;
        else
            [w(1,face),grads(face,:)] = mean_value_2D(verts(face,:), xin);
            varargout{1} = w;
            varargout{2} = grads;
        end
        return
    end
end
for f=1:length(faces)
    face = faces{f};
    fc = mean(verts(face,:));
    uc = (fc - xin) / norm(fc - xin);
    for ff=1:length(face)
        if ff == length(face)
            v = [face(ff),face(1)];
        else
            v = [face(ff),face(ff+1)];
        end
        l1 = norm(u(v(2),:) - uc);
        l2 = norm(uc - u(v(1),:));
        l3 = norm(u(v(2),:) - u(v(1),:));
        t1 = 2*asin(l1/2);
        t2 = 2*asin(l2/2);
        t3 = 2*asin(l3/2);
        h = sum([t1,t2,t3])/2;
        c1 = 2*sin(h)*sin(h-t1) / (sin(t2)*sin(t3)) - 1;
        c2 = 2*sin(h)*sin(h-t2) / (sin(t3)*sin(t1)) - 1;
        c3 = 2*sin(h)*sin(h-t3) / (sin(t2)*sin(t1)) - 1;
        s1 = sign(det([u(v,:)',uc'])) * sqrt(1-c1^2);
        s2 = sign(det([u(v,:)',uc'])) * sqrt(1-c2^2);
        s3 = sign(det([u(v,:)',uc'])) * sqrt(1-c3^2);
        if abs(s1) < 1e-14 || abs(s2) < 1e-14 || abs(s3) < 1e-14
            % point lies outside the spherical shell of the projected triangle
            continue
        end
        w(1,v(1)) = w(1,v(1)) + (t1 - c2*t3 - c3*t2) / (d(v(1))*sin(t2)*s3);
        w(1,v(2)) = w(1,v(2)) + (t2 - c3*t1 - c1*t3) / (d(v(2))*sin(t3)*s1);
    end
end
w = w / sum(w);
varargout{1} = w;
if nout == 2, varargout{2} = grads; end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_planes(verts, faces)
nf = length(faces);
out = cell(nf, 1);
for i=1:nf
    f = faces{i};
    out{i} = createPlane(verts(f,:));
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

