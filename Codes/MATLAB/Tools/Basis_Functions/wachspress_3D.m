function varargout = wachspress_3D(verts, xin, faces, varargin)
% Get input/output information
% ----------------------------
if nargin <= 1, error('No inputs specified.'); end
nv = size(verts,1);
nx = size(xin,1);
dim1 = size(verts,2);
dim2 = size(xin,2);
if dim1 ~= dim2, error('Dimensions between inputs does not match.'); end

% Get Wachpress interpolant values/grads
% --------------------------------------
p = get_face_planes(verts, faces);
if nargout == 1
    varargout{1} = get_values(verts, faces, xin, nv, nx);
elseif nargout == 2
    [varargout{1},varargout{2}] = get_values(verts, faces, xin, nv, nx, dim1);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Function Lists
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = get_values(verts, faces, xin, nv, nx, dim)
nout = nargout;
vals = zeros(nx, nv);
[g,un,f2face] = get_facets(verts, faces);
if nout == 2
    grads = zeros(nx, nv, dim);
end
if nout == 1
    for j=1:nx
        vals(j,:) = get_single_value(verts, faces, xin(j,:), g, un, f2face, nv);
    end
else
    for j=1:nx
        [vals(j,:), grads(j,:,:)] = get_single_value(verts, faces, xin(j,:), g, un, f2face, nv);
    end
end
% outputs
varargout{1} = vals;
if nout == 2, varargout{2} = grads; end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = get_single_value(verts, faces, xin, g, un, f2face, nv)
nout = nargout;
if nout == 2
    R = zeros(nv,3);
    grads = zeros(nv,3);
end
w = zeros(1,nv);
xx = verts - ones(nv,1)*xin;
for i=1:nv
    if norm(xx(i,:)) < 1e-14
        w(i) = 1;
        varargout{1} = w;
        if nout==2
            grads = nan;
            varargout{2} = grads; 
        end
        return
    end
end
for i=1:nv
    f = g{i};
    k = length(f);
    p = zeros(k,3);
    for j=1:k
        h = dot(verts(i,:) - xin, un(f(j),:));
        if abs(h) < 1e-14
            face = f2face(f(j));
            if nout == 1
                w(1,faces{face}) = wachpress_2D(verts(faces{face},:), xin);
                varargout{1} = w;
            else
                [w(1,faces{face}), grads(faces{face},:)] = wachpress_2D(verts(faces{face},:), xin);
                varargout{1} = w;
                varargout{2} = grads;
            end
            return
        end
        p(j,:) = un(f(j),:) / h;
    end
    wloc = zeros(k-2,1);
    if nout==2, Rloc = zeros(k-2,3); end
    for j=1:k-2
        wloc(j) = det([p(j,:); p(j+1,:); p(k,:)]);
        if nout==2, Rloc(j,:) = p(j,:) + p(j+1,:) + p(k,:); end
    end
    w(i) = sum(wloc);
    if nout==2, R(i,:) = (wloc' * Rloc) / w(i); end
end
varargout{1} = w / sum(w);
if nout == 2
    phiR = varargout{1}' * R;
    for d=1:3
        varargout{2}(:,d) = varargout{1} .* (R(:,d) - phiR(:,d));
    end
    varargout{2} = grads;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,f2face] = get_vertex_connections(verts, faces, nv)
xc = mean(verts);
out = cell(nv,1);
count = 0;
for f=1:length(faces)
    face = faces{f};
    for j=1:length(face)
        count = count + 1;
        if j==1
            v = [length(face),j+1];
        elseif j==length(face)
            v = [j-1,1];
        else
            v = [j-1,j+1];
        end
        vert = face(j);
        out{vert} = [out{vert},face(v)];
        f2face(count) = f;
    end
end
for i=1:nv
    out{i} = unique(out{i});
end
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
function [graphs, facet_norms, f2face] = get_facets(verts, faces)
vc = mean(verts);
nf = length(faces);
nv = size(verts,1);
graphs = cell(nv,1);
[gg,f2face] = get_vertex_connections(verts, faces, nv);
tf = 0;
for i=1:nf
    f = faces{i};
    tf = tf + length(f);
end
facet_norms = zeros(tf,3);
count = 0;
for i=1:nv
    vf = gg{i};
    for j=1:length(vf)
        count = count + 1;
        if j==length(vf)
            e = [j,1];
        else
            e = [j,j+1];
        end
        vvf = vf(e);
        tvecs = verts(vvf,:) - ones(size(verts(vvf,:),1),1)*verts(i,:);
        facet_norms(count,:) = cross(tvecs(2,:),tvecs(1,:));
        if dot(verts(i,:)-vc, facet_norms(count,:)) < 0
            facet_norms(count,:) = -1*facet_norms(count,:);
        end
        graphs{i} = [graphs{i},count];
    end
end
return