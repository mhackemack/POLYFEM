%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D Lagrange Basis Function Generator - Volume Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the Triangular and Quadrilateral basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = Lagrange_volume(varargin)
nin = nargin;
nout = nargout;
%
% Collect Input Arguments
% -----------------------
nverts = varargin{5};
verts = varargin{1};
flags = varargin{3};
deg = varargin{4};
%
% Perform Error Checking Operations
% ---------------------------------
[m,n] = size(verts);
if n>m, verts = verts'; end
[nv,dim] = size(verts);
% if nargin < 4, flags = [1,1,1]; end
if nout ~= sum(flags), error('Number of outputs does not match input flags.'); end
if dim == 1
    if nverts > 2, error('Only 2 points in 1D.'); end
    ctype = 1;
elseif dim == 2
    if nverts ~= 3 && nverts ~= 4, error(''); end
    if nverts == 3, ctype = 1; end
    if nverts == 4, ctype = 2; end
elseif dim == 3
    if nverts ~= 4 && nverts ~= 8, error(''); end
    if nverts == 4, ctype = 1; end
    if nverts == 8, ctype = 2; end
end
%
% Get Jacobian and Other Information
% ----------------------------------
nbf = get_num_basis_funcs(ctype, dim, deg);
if flags(2) || flags(3) || ctype == 2
    [qx, qw] = get_quadrature(ctype, dim);
    nqx = length(qw);
    db = evaluate_ref_gradients(ctype, dim, deg, qx, nqx, nbf);
end
if flags(3) || ctype == 2
    b = evaluate_ref_vals(ctype, dim, deg, qx);
end
if ctype == 1
    [~, invJ, detJ, vol] = get_jacobian(verts, dim);
else
    J = evaluate_numerical_jacobian(dim, verts, qx, deg);
    detJ = zeros(nqx,1);
    for q=1:nqx
        detJ(q) = det(J{q});
    end
    if flags(2) || flags(3)
        invJ = cell(nqx,1);
        for q=1:nqx
            invJ{q} = inv(J{q});
        end
    end
end
% 
% Switch Cell Type and Degree
% ---------------------------
counter = 1;
% mass matrix
if flags(1)
    if ctype == 1
        varargout{counter} = vol*get_ref_mass_matrix(dim, deg);
        counter = counter + 1;
    else
        M = zeros(nbf, nbf);
        for q=1:nqx
            M = M + (b(q,:)'*detJ(q)*b(q,:))*qw(q);
        end
        varargout{counter} = M;
        counter = counter + 1;
    end
end
switch(ctype)
    case(1)
        % stiffness matrix
        if flags(2)
            K = zeros(nbf, nbf);
            tJ = (invJ*invJ')*detJ;
            for q=1:nqx
                K = K + (db(:,:,q)*tJ*db(:,:,q)')*qw(q);
            end
            varargout{counter} = K;
            counter = counter + 1;
        end
        % gradient matrix
        if flags(3)
            tJ = invJ*detJ;
            G = cell(dim, 1);
            for d=1:dim 
                G{d} = zeros(nbf, nbf);
            end
            for q=1:nqx
                tdb = db(:,:,q)*tJ;
                for d=1:dim
                    G{d} = G{d} + (tdb(:,d)*b(q,:))*qw(q);
                end
            end
            varargout{counter} = G;
        end
    case(2)
        % stiffness matrix
        if flags(2)
            K = zeros(nbf, nbf);
            for q=1:nqx
                tJ = (invJ{q} * invJ{q}') * detJ(q);
                K = K + (db(:,:,q)*tJ*db(:,:,q)')*qw(q);
            end
            varargout{counter} = K;
            counter = counter + 1;
        end
        % gradient matrix
        if flags(3)
            G = cell(dim, 1);
            for d=1:dim 
                G{d} = zeros(nbf, nbf);
            end
            for q=1:nqx
                tdb = db(:,:,q)*invJ{q}*detJ(q);
                for d=1:dim
                    G{d} = G{d} + (tdb(:,d)*b(q,:))*qw(q);
                end
            end
            varargout{counter} = G;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, invJ, detJ, vol] = get_jacobian(verts, dim)
J = zeros(dim);
vverts = verts';
for d=1:dim
    J(:,d) = vverts(:,d+1) - vverts(:,1);
end
invJ = inv(J);
detJ = det(J);
if dim == 1
    vol = detJ;
else
    vol = detJ / (dim * (dim - 1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = get_ref_mass_matrix(dim, deg)
if dim == 1
    if deg == 1
        m = [2,1;1,2]./6;
    elseif deg == 2
        m = [4,-1,2;-1,4,2;2,2,16]./30;
    elseif deg == 3
        m = [8/105,19/1680,33/560,-3/140;19/1680,8/105,-3/140,33/560;33/560,-3/140,27/70,-27/560;-3/140,33/560,-27/560,27/70];
    end
else
    if deg == 1
        if dim==2
            m = [2,1,1;1,2,1;1,1,2]./12;
        elseif dim==3
            m = [2,1,1,1;1,2,1,1;1,1,2,1;1,1,1,2]./20;
        end
    elseif deg == 2
        if dim==2
            m = [6,-1,-1,0,-4,0;...
                 -1,6,-1,0,0,-4;...
                 -1,-1,6,-4,0,0;...
                 0,0,-4,32,16,16;...
                 -4,0,0,16,32,16;...
                 0,-4,0,16,16,32]./180;
        elseif dim==3
            
        end
    elseif deg == 3
        if dim == 2
            
        elseif dim == 3
            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_ref_vals(ctype, dim, deg, qx)
fun = get_basis_funcs(ctype, dim);
out = fun(deg, qx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_ref_gradients(ctype, dim, deg, qx, nqx, n)
out = zeros(n,dim,nqx);
fun = get_basis_grad_funcs(ctype, dim);
for q=1:nqx
    out(:,:,q) = fun(deg, qx(q,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_num_basis_funcs(ctype, dim, deg)
if dim == 1
    out = 1 + deg;
elseif dim == 2
    if ctype == 1
        if deg == 3
            out = 10;
        else
            out = 3*deg;
        end
    else
        out = (deg+1)^2;
    end
else
    if ctype == 1
        if deg == 1
            out = 4;
        elseif deg == 2
            out = 10;
        end
    else
        out = (deg+1)^3;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_funcs(ctype, dim)
out = get_lagrange_function('vals', ctype, dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_grad_funcs(ctype, dim)
out = get_lagrange_function('grads', ctype, dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_numerical_jacobian(dim, verts, x, deg)
n = size(x,1);
out = cell(n,1);
% Evaluate Jacobian
if dim == 2
    db = eval_quad_grads_for_jac(1, x);
    for i=1:n
        out{i} = [db(:,1,i)'*verts(:,1), db(:,2,i)'*verts(:,1);...
                  db(:,1,i)'*verts(:,2), db(:,2,i)'*verts(:,2)];
    end
elseif dim == 3
    db = eval_hex_grads_for_jac(1, x);
    for i=1:n
        out{i} = [db(:,1,i)'*verts(:,1), db(:,2,i)'*verts(:,1), db(:,3,i)'*verts(:,1);...
                  db(:,1,i)'*verts(:,2), db(:,2,i)'*verts(:,2), db(:,3,i)'*verts(:,2);...
                  db(:,1,i)'*verts(:,3), db(:,2,i)'*verts(:,3), db(:,3,i)'*verts(:,3)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_quadrature(ctype, dim)
if ctype == 1
    if dim == 1
        [qx, qw] = get_legendre_gauss_quad(4);
    elseif dim == 2
        [qx, qw] = Quad_On_Triangle(6);
    elseif dim == 3
        [qx, qw] = Quad_On_Tetra(10);
    end
else
    [tqx, tqw] = get_legendre_gauss_quad(4); nt = length(tqw);
    if dim == 2
        qw = tqw*tqw'; qw = qw(:); qw = qw / sum(qw);
        qxx = zeros(nt, nt); qxy = zeros(nt, nt);
        for i=1:nt
            qxx(i,:) = tqx(i);
            qxy(:,i) = tqx(i);
        end
        qx = [qxx(:), qxy(:)];
    else
        qxx = zeros(nt, nt, nt); qxy = zeros(nt, nt, nt); qxz = zeros(nt, nt, nt); 
        qw = zeros(nt, nt, nt); qqw = tqw*tqw';
        for i=1:nt
            qxx(i,:,:) = tqx(i);
            qxy(:,i,:) = tqx(i);
            qxz(:,:,i) = tqx(i);
            qw(:,:,i)  = qqw*tqw(i);
        end
        qx = [qxx(:), qxy(:), qxz(:)];
        qw = qw(:); qw = qw / sum(qw);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_numbering(ctype, dim)
if dim == 2
    if ctype == 1
        out = cell(3,1);
        out{1} = [1,2];
        out{2} = [2,3];
        out{3} = [3,1];
    else
        out = cell(4,1);
        out{1} = [1,2];
        out{2} = [2,3];
        out{3} = [3,4];
        out{4} = [4,1];
    end
elseif dim == 3
    if ctype == 1
        out = cell(4,1);
        out{1} = [1,2,3,4];
    else
        out = cell(6,1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_quad_hex_volume(verts, ctype, dim)
out = 0;
vcenter = mean(verts);
faces = get_face_numbering(ctype, dim);
for i=1:length(faces)
    if dim == 2
        out = out + heron_triangle_area([verts(faces{i},:);vcenter]);
    else
        fcenter = mean(verts(faces{i},:));
        for j=1:length(faces{i})
            out = out + tet_volume([verts(faces{i},:);fcenter;vcenter]);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function area = heron_triangle_area(v)
a = norm(v(2,:) - v(1,:));
b = norm(v(3,:) - v(2,:));
c = norm(v(1,:) - v(3,:));
s = (a+b+c)/2;
area = sqrt(s*(s-a)*(s-b)*(s-c));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vol = tet_volume(verts)
J = zeros(3,3);
vverts = verts';
for d=1:dim
    J(:,d) = vverts(:,d+1) - vverts(:,1);
end
detJ = det(J);
vol = detJ / 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = eval_hex_grads_for_jac(deg, xx)
n = size(xx,1); out = zeros((deg+1)^3,3,n);
for i=1:n
x = xx(i,1); y = xx(i,2); z = xx(i,3);
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
function out = eval_quad_grads_for_jac(deg, xx)
n = size(xx,1); out = zeros((deg+1)^2,2,n);
for i=1:n
    s = xx(i,1); t = xx(i,2);
    if deg == 1
        out(:,:,i) = [t-1,s-1;1-t,-s;t,s;-t,1-s];
    elseif deg == 2
        out(:,:,i) = [  8*s*t*t-12*s*t+4*s-6*t*t+9*t-3,    8*s*s*t-6*s*s-12*s*t+9*s+4*t-3;...
                 8*s*t*t-12*s*t+4*s-2*t*t+3*t-1,    8*s*s*t-6*s*s-4*s*t+3*s;...
                 8*s*t*t-4*s*t-2*t*t+t,             8*s*s*t-2*s*s-4*s*t+s;...
                 8*s*t*t-4*s*t-6*t*t+3*t,           8*s*s*t-2*s*s-12*s*t+3*s+4*t-1;...
               -16*s*t*t+24*s*t-8*s+8*t*t-12*t+4, -16*s*s*t+12*s*s+16*s*t-12*s;...
               -16*s*t*t+16*s*t+4*t*t-4*t,        -16*s*s*t+8*s*s+8*s*t-4*s;...
               -16*s*t*t+8*s*t+8*t*t-4*t,         -16*s*s*t+4*s*s+16*s*t-4*s;...
               -16*s*t*t+16*s*t+12*t*t-12*t,      -16*s*s*t+8*s*s+24*s*t-12*s-8*t+4;...
                32*s*t*t-32*s*t-16*t*t+16*t,       32*s*s*t-16*s*s-32*s*t+16*s];
    elseif deg == 3

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%