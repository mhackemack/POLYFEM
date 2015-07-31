%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Lagrange Quadrature Generator for Simplexes and
%                   Quads/Hexes
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
function varargout = Lagrange_quad_gen(varargin)
nin = nargin;
nout = nargout;
%
% Collect Input Arguments
% -----------------------
nverts = varargin{5};
verts = varargin{1};
q_ord = varargin{2};
deg = varargin{4};
%
% Perform Error Checking Operations
% ---------------------------------
[m,n] = size(verts);
if n>m, verts = verts'; end
[nv,dim] = size(verts);
if dim == 1
    ctype = 1;
elseif dim == 2
    if nv == 3
        ctype = 1;
    elseif nv == 4
        ctype = 2;
    end
else
    if nv == 4
        ctype = 1;
    elseif nv == 8
        ctype = 2;
    end
end
% Output argument checking
if nargout < 2
    error('Too few output arguments - need at least 2.')
elseif nargout < 3
    mass_out = 0;
    grad_out = 0;
elseif nargout == 3
    mass_out = 1;
    grad_out = 0;
elseif nargout == 4
    mass_out = 1;
    grad_out = 1;
end
%
% Get Jacobian and Other Information
% ----------------------------------
nbf = get_num_basis_funcs(ctype, dim, deg);
[qx, qw] = get_quadrature(ctype, dim, q_ord); nqx = length(qw);
[J, invJ, detJ] = get_jacobian(ctype, verts, dim, qx);
%
% Get Reference Values and Allocate Memory
% ----------------------------------------
quad_nodes_out = zeros(nqx, dim);
quad_wts_out = qw.*detJ;
if mass_out
    basis_vals_out = evaluate_ref_vals(ctype, dim, deg, qx);
end
if grad_out
    basis_grads_out = evaluate_ref_gradients(ctype, dim, deg, qx, nqx, nbf);
end
% Set global quadrature node locations
for q=1:nqx
    quad_nodes_out(q,:) = verts(1,:) + (J{q}*qx(q,:)')';
end
% Set Outputs
varargout{1} = quad_nodes_out;
varargout{2} = quad_wts_out;
if mass_out
    varargout{3} = basis_vals_out;
end
if grad_out
    if mass_out
        varargout{4} = basis_grads_out;
    else
        varargout{3} = basis_grads_out;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, invJ, detJ] = get_jacobian(ctype, verts, dim, qx)
n = size(qx, 1);
J = cell(n,1); invJ = cell(n,1); detJ = zeros(n,1);
if ctype == 1
    JJ = zeros(dim);
    vverts = verts';
    for d=1:dim
        JJ(:,d) = vverts(:,d+1) - vverts(:,1);
    end
    dJ = det(JJ);
    iJ = inv(JJ);
    for i=1:n
        J{i} = JJ;
        detJ(i) = dJ;
        invJ{i} = iJ;
    end
elseif ctype == 2
    J = evaluate_numerical_jacobian(dim, verts, qx);
    for i=1:n
        detJ(i) = det(J{i});
        invJ{i} = inv(J{i});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_numerical_jacobian(dim, verts, x)
n = size(x,1);
out = cell(n,1);
% Evaluate Jacobian
if dim == 2
    xall = verts(1,1) - verts(2,1) + verts(3,1) - verts(4,1);
    yall = verts(1,2) - verts(2,2) + verts(3,2) - verts(4,2);
    x21 = verts(2,1) - verts(1,1); y21 = verts(2,2) - verts(1,2);
    x41 = verts(4,1) - verts(1,1); y41 = verts(4,2) - verts(1,2);
    for i=1:n
        out{i} = [x21 + xall*x(i,2), x41 + xall*x(i,1);y21 + yall*x(i,2), y41 + yall*x(i,1)];
    end
elseif dim == 3
    db = eval_hex_grads_for_jac(x);
    for i=1:n
        out{i} = [db(:,1,i)'*verts(:,1), db(:,2,i)'*verts(:,1), db(:,3,i)'*verts(:,1);...
                  db(:,1,i)'*verts(:,2), db(:,2,i)'*verts(:,2), db(:,3,i)'*verts(:,2);...
                  db(:,1,i)'*verts(:,3), db(:,2,i)'*verts(:,3), db(:,3,i)'*verts(:,3)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_quadrature(ctype, dim, q_ord)
if ctype == 1
    if dim == 1
        [qx, qw] = get_legendre_gauss_quad(max(q_ord,6));
    elseif dim == 2
        [qx, qw] = Quad_On_Triangle(max(q_ord,6));
    elseif dim == 3
        [qx, qw] = Quad_On_Tetra(max(q_ord,10));
    end
else
    [tqx, tqw] = get_legendre_gauss_quad(max(q_ord,6)); nt = length(tqw);
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
            qxy(:,:,i) = tqx(i);
            qw(:,:,i)  = qqw*tqw(i);
        end
        qx = [qxx(:), qxy(:), qxz(:)];
        qw = qw(:); qw = qw / sum(qw);
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
function out = get_basis_funcs(ctype, dim)
out = get_lagrange_function('vals', ctype, dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_grad_funcs(ctype, dim)
out = get_lagrange_function('grads', ctype, dim);
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
function out = get_num_basis_funcs(ctype, dim, deg)
if dim == 1
    out = 2 + (deg-1);
elseif dim == 2
    if ctype == 1
        out = 3*deg;
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