%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          General Quadrature Generator - Volume
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce a quadrature set for general
%                   polygons/polyhedra.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_general_volume_quadrature(verts, faces, deg, poly_bool)
% Get General Input Information 
% ------------------------------------------------------------------------------
[nverts, dim] = size(verts);
if dim == 1
    [qx, qw] = lgwt(deg, verts(1), verts(2));
    return
elseif dim == 2
    nfaces = 1; faces{1} = 1:nverts;
elseif dim == 3
    nfaces = length(faces);
end
rcenter = mean(verts);
ctype = get_cell_type(nverts, dim, nfaces);
nvs = get_number_side_volumes(dim, ctype, faces, poly_bool);
% Get Quadrature Prelims
% ------------------------------------------------------------------------------
[rqx, rqw] = get_reference_quadrature(ctype, dim, deg, poly_bool);
nrqx = length(rqw); nqx  = nvs*nrqx;
qx   = zeros(nqx, dim);
qw   = zeros(nqx,1);
% Get Quadrature Set
% ------------------------------------------------------------------------------
if ctype == 3 || poly_bool
    for f=1:nfaces
        ff = faces{f}; nf = length(ff);
        fcenter = mean(verts(ff,:));
        cind = 0;
        for i=1:nf
            if i==nf
                ii = [i,1];
            else
                ii = [i,i+1];
            end
            ffi = ff(ii);
            if dim == 2
                vv = [verts(ffi,:);rcenter];
            elseif dim == 3
                vv = [verts(ffi,:);fcenter;rcenter];
            end
            J = get_simplex_jacobian(dim,vv);
            detJ = det(J);
            svol = detJ / (dim * (dim-1));
            v0 = vv(1,:);
            for q=1:nrqx
                cind = cind + 1;
                qx(cind,:) = v0 + (J*rqx(q,:)')';
                qw(cind) = rqw(q) * svol;
            end
        end
    end
elseif ctype == 1
    J = get_simplex_jacobian(dim, verts);
    svol = det(J) / factorial(dim);
    v0 = verts(1,:);
    for q=1:nrqx
        qx(q,:) = v0 + (J*rqx(q,:)')';
    end
    qw = rqw * svol;
else
    v0 = verts(1,:);
    [J, detJ] = get_quad_hex_jacobians(verts, rqx);
    for q=1:nrqx
        qx(q,:) = v0 + (J{q}*rqx(q,:)')';
    end
    qw = rqw.*detJ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctype = get_cell_type(nverts, dim, nfaces)
if dim == 1
    ctype = 1;
elseif dim == 2
    if nverts == 3
        ctype = 1;
    elseif nverts == 4
        ctype = 2;
    else
        ctype = 3;
    end
elseif dim == 3
    if nfaces == 4
        ctype = 1;
    elseif nfaces == 6
        ctype = 2;
    else
        ctype = 3;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_number_side_volumes(dim, ctype, faces, pb)
if ctype == 3 || pb
    if dim == 2
        out = length(faces{1});
    else
        nf = length(faces);
        out = 0;
        for f=1:nf
            out = out + length(faces{f});
        end
    end
else
    out = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_reference_quadrature(ctype, dim, deg, pb)
if pb || ctype ~= 2
    if dim == 2
        [qx, qw] = Quad_On_Triangle(deg);
    else
        [qx, qw] = Quad_On_Tetra(deg);
    end
    qw = qw / sum(qw);
    return
else
    [qx, qw] = get_cart_ref_quad(dim, deg);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_cart_ref_quad(dim, deg)
[tqx, tqw] = get_legendre_gauss_quad(deg); nt = length(tqw);
% Switch between dimension
switch(dim)
    case(1)
        qx = tqx; qw = tqw;
    case(2)
        qw = tqw*tqw'; qw = qw(:); qw = qw / sum(qw);
        qxx = zeros(nt, nt); qxy = zeros(nt, nt);
        for i=1:nt
            qxx(i,:) = tqx(i); qxy(:,i) = tqx(i);
        end
        qx = [qxx(:), qxy(:)];
    case(3)
        qxx = zeros(nt, nt, nt); qxy = zeros(nt, nt, nt); qxz = zeros(nt, nt, nt); 
        qw = zeros(nt, nt, nt); qqw = tqw*tqw';
        for i=1:nt
            qxx(i,:,:) = tqx(i); qxy(:,i,:) = tqx(i); qxz(:,:,i) = tqx(i);
            qw(:,:,i)  = qqw*tqw(i);
        end
        qx = [qxx(:), qxy(:), qxz(:)];
        qw = qw(:); qw = qw / sum(qw);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, detJ] = get_quad_hex_jacobians(verts, x)
[n, dim] = size(x);
J = cell(n, 1);
detJ = zeros(n, 1);
if dim == 2
    xall = verts(1,1) - verts(2,1) + verts(3,1) - verts(4,1);
    yall = verts(1,2) - verts(2,2) + verts(3,2) - verts(4,2);
    x21 = verts(2,1) - verts(1,1); y21 = verts(2,2) - verts(1,2);
    x41 = verts(4,1) - verts(1,1); y41 = verts(4,2) - verts(1,2);
    for i=1:n
        J{i} = [x21 + xall*x(i,2), x41 + xall*x(i,1);y21 + yall*x(i,2), y41 + yall*x(i,1)];
        JJ = J{i};
        detJ(i) = JJ(1,1)*JJ(2,2)-JJ(2,1)*JJ(1,2);
    end
elseif dim == 3
    db = eval_hex_grads_for_jac(x);
    for i=1:n
        J{i} = [db(:,1,i)'*verts(:,1), db(:,2,i)'*verts(:,1), db(:,3,i)'*verts(:,1);...
                db(:,1,i)'*verts(:,2), db(:,2,i)'*verts(:,2), db(:,3,i)'*verts(:,2);...
                db(:,1,i)'*verts(:,3), db(:,2,i)'*verts(:,3), db(:,3,i)'*verts(:,3)];
        JJ = J{i};
        detJ(i) = JJ(1,1)*(JJ(3,3)*JJ(2,2)-JJ(3,2)*JJ(2,3)) - JJ(2,1)*(JJ(3,3)*JJ(1,2)-JJ(3,2)*JJ(1,3)) + JJ(3,1)*(JJ(2,3)*JJ(1,2)-JJ(2,2)*JJ(1,3));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = eval_hex_grads_for_jac(xx)
x = xx(:,1); y = xx(:,2); z = xx(:,3);
n = size(xx,1); out = zeros(8,3,n);
for i=1:n
out(:,:,i) = [  -y.*z+y+z-1,  -x.*z+x+z-1,  -x.*y+x+y-1;...
                (1-y).*(1-z),  x.*z-x,       x.*y-x;...
                y.*(1-z),      x.*(1-z),    -x.*y;...
                y.*z-y,       (1-x).*(1-z),  x.*y-y;...
                y.*z-z,       x.*z-z         (1-x).*(1-y);...
                (1-y).*z,     -x.*z          x.*(1-y);...
                y.*z,          x.*z,         x.*y;...
                -y.*z,         (1-x).*z,     (1-x).*y];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%