function varargout = convex_poly_quad( v, levels, faces)
% input/output
nout = nargout;
if size(v,1) < size(v,2)
    v = v';
end
[nv, dim] = size(v);
% error checking
if nout ~= 2 && nout ~= 4
    error('Incorrect output structure.')
end
if dim == 3 && nout == 4
    if nargin < 3
        error('Missing Face Structure.')
    end
end

% convex hull volume information
[K, V] = convhulln(v);

switch(dim)
    case(2)
        if nv == 3
            [vx, vw] = triquad(levels,v);
        elseif nv == 4
            [vx,vw] = get_unit_square_quad(levels);
            vx = map_hypercube_to_hyperrect(dim,v,vx);
        else
            [vx, vw] = get_2d_quad(v, levels);
        end
        if nout == 4
            fx = cell(nv, 1);
            fw = cell(nv, 1);
            for i=1:nv
                if i==nv
                    vv = [i,1];
                else
                    vv = [i,i+1];
                end
                [x, w] = get_1d_quad(v(vv,:), levels);
                fx{i} = [fx{i};x];
                fw{i} = [fw{i};w];
            end
        end
    case(3)
        if nv == 4
            [vx, vw] = tetraquad(levels,v);
        else 
            [vx, vw] = get_3d_quad(v, K, levels);
        end
        if nout == 4
            nf = length(faces);
            fx = cell(nf,1);
            fw = cell(nf,1);
            for f=1:nf
                [x, w] = get_2d_quad(v(faces{f},:), levels);
                fx{f} = [fx{f};x];
                fw{f} = [fw{f};w];
            end
        end
end

% Set outputs
vw = vw * V / sum(vw);
varargout{1} = vx;
varargout{2} = vw;
if nout == 4
    varargout{3} = fx;
    varargout{4} = fw;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = get_1d_quad(v, levels)
dist = norm(v(2,:) - v(1,:));
[xx,w] = lgwt(levels, 0, dist);
x = v(1,:)*ones(size(xx,1),1) + xx;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx,ww] = get_2d_quad(verts, levels)
vc = mean(verts);
[nv, dim] = size(verts);
qx = []; ww = [];
for i=1:nv
    if i==nv
        v = [i,1];
    else
        v = [i,i+1];
    end
    tverts = flipud([verts(v,:);vc]);
    if dim == 2
        [x, w] = triquad(levels, tverts);
        qx = [qx;x];
        ww = [ww;w];
    else
        tform = transform_3d_plane_to_2d(tverts);
        
        J = get_jacobian(tverts);
        qx = [qx;x];
        ww = [ww;w];
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx,ww] = get_3d_quad(v, tets, levels)
vc = mean(v);
nv = size(v,1);
nt = size(tets, 1);
qx = []; ww = [];
for i=1:nt
    vv = [v(tets(i,:),:);vc];
    [x, w] = tetraquad(levels, vv);
    qx = [qx;x];
    ww = [ww;w];
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = get_jacobian(v)
[nv, dim] = size(v);
J = zeros(dim,dim);
for i=1:dim
    J(:,i) = v(i+1,:)' - v(1,:)';
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = get_unit_square_quad(levels)
[x,w]=lgwt(levels, 0 ,1);
[x,y]=meshgrid(x,x);
x=[x(:),y(:)];
w = w*w'; w=w(:);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,wout] = get_unit_cube_quad(levels)
[x,w]=lgwt(levels, 0 ,1);
[x,y,z]=meshgrid(x,x,x);
x=[x(:),y(:),z(:)];
ww = w*w'; www=ww(:); wout = [];
for i=1:levels
    tw = www*w(i);
    wout = [wout;tw];
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = map_hypercube_to_hyperrect(dim, v, xin)
A = get_linear_map(dim);
out = zeros(size(xin,1),dim);
if dim == 2
    a = A\v(:,1);
    b = A\v(:,2);
    for i=1:size(xin,1)
        out(i,1) = a(1) + a(2)*xin(i,1) + a(3)*xin(i,2) + a(4)*xin(i,1)*xin(i,2);
        out(i,2) = b(1) + b(2)*xin(i,1) + b(3)*xin(i,2) + b(4)*xin(i,1)*xin(i,2);
    end
elseif dim == 3
    
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_linear_map(dim)
if dim == 2
    out = [1,0,0,0;1,1,0,0;1,1,1,1;1,0,1,0];
elseif dim == 3
    
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = map_3d_triangle_to_2d(xin)
out = zeros(3,2);
l1 = norm(xin(2,:) - xin(1,:));
l2 = norm(xin(3,:) - xin(2,:));
l3 = norm(xin(1,:) - xin(3,:));
a1 = acos(dot(xin(3,:)-xin(1,:),xin(2,:)-xin(1,:)) / (l1*l3));
a2 = acos(dot(xin(3,:)-xin(2,:),xin(1,:)-xin(2,:)) / (l1*l2));
a3 = acos(dot(xin(1,:)-xin(3,:),xin(2,:)-xin(3,:)) / (l2*l3));

out(2,:) = [l1,0];
out(3,:)
return