%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Cartesian Lagrange Main Generation Function
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB function to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the Lagrange basis functions on quad and hex 
%                   meshes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input Space:    1) Number of geometric vertices
%                   2) Vertices
%                   3) Face Vertices
%                   4) FEM Order
%                   5) FEM Lumping Boolean
%                   6) Volumetric Matrix Flags
%                   7) Surface Matrix Flags
%                   8) Quadrature boolean
%                   9) Quadrature Order (Optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = bf_cell_func_Lagrange_Cartesian( varargin )
% Collect Input/Output Arguments
% ------------------------------------------------------------------------------
nout = nargout;
nverts = varargin{1};
verts = varargin{2}(1:nverts,:);
faces = varargin{3}; nf = length(faces);
ord = varargin{4};
lump_bool = varargin{5};
v_flags = varargin{6};
s_flags = varargin{7};
q_bool = varargin{8};
q_ord = ord+2;
if nargin > 8
    if ~isempty(varargin{9}),q_ord = varargin{9};end
end
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
% Compute and exit immediately if 1D
% ------------------------------------------------------------------------------
if dim == 1
    [bf_V,bf_S,QV,QS] = bf_func_1D(varargin{:});
    varargout = {bf_V, bf_S, QV, QS};
    return
end
% ------------------------------------------------------------------------------
% Quick Error Checking
% ------------------------------------------------------------------------------
if dim == 2
    if ord > 3, error('Cannot go higher than 3rd order in 3D at this time.'); end
    if nv > 4, error('Polygons are not allowed with this basis function.'); end
else
    if ord > 1, error('Cannot go higher than 1st order in 3D at this time.'); end
end
% Generate Basis Function Space
% ------------------------------------------------------------------------------
% Retrieve Reference Function Handles
% -----------------------------------
ref_val_func = get_lagrange_function('vals', 2, dim);
ref_grad_func = get_lagrange_function('grads', 2, dim);
% Generate Reference Quadrature Spaces
% ------------------------------------
ntot = get_num_dofs(dim, ord); zt = ones(1,ntot);
f_dofs = get_face_dofs(dim, ord);
if dim == 2
    ref_x = [0,0;1,0;1,1;0,1];
    svol = polygonArea(verts);
elseif dim == 3
    ref_x = [0,0,0;1,0,0;1,1,0;0,1,0;0,0,1;1,0,1;1,1,1;0,1,1];
    [~,svol] = convhull(verts);
end
[rqx_V, rqw_V] = get_ref_quadrature(dim, q_ord);   nq_V = length(rqw_V);
[rqx_S, rqw_S] = get_ref_quadrature(dim-1, q_ord); nq_S = length(rqw_S);
[trqx_S, trqw_S] = get_surface_ref_quadrature(dim, rqx_S, rqw_S, ref_x);
zs = ones(nq_S, 1);
% Compute Jacobian
% ----------------
[J,invJ] = evaluate_numerical_jacobian(dim, verts, rqx_V);
JS = cell(nf,1); invJS = cell(nf,1);
for f=1:nf
    [JS{f},invJS{f}] = evaluate_numerical_jacobian(dim, verts, trqx_S{f});
end
% Compute Real Space Quadrature and Basis Values/Gradients
% --------------------------------------------------------
% Volume
qx_V = zeros(nq_V, dim); qw_V = rqw_V*svol;
bvals_v = ref_val_func(ord, rqx_V);
bgrads_v = zeros(ntot,dim,nq_V);
x0 = verts(1,:);
for q=1:nq_V
    qx_V(q,:) = x0 + (J{q}*rqx_V(q,:)')';
    tg = ref_grad_func(ord, rqx_V(q,:));
    bgrads_v(:,:,q) = tg*invJ{q};
end
% Surface
bvals_s = cell(nf, 1); bgrads_s = cell(nf, 1);
qx_S    = cell(nf, 1); qw_S     = cell(nf, 1);
for f=1:nf
    bgrads_s{f} = zeros(ntot, dim, nq_S);
    fv = faces{f};
    fverts = verts(fv,:);
    f0 = fverts(1,:);
    if dim == 2
        dx = diff(fverts);
        len = norm(dx);
        qw_S{f} = rqw_S*len;
        qx_S{f} = zs*f0 + rqx_S*dx;
        bvals_s{f} = ref_val_func(ord,trqx_S{f});
        for q=1:nq_S
            bgrads_s{f}(:,:,q) = ref_grad_func(ord, trqx_S{f}(q,:))*invJS{f}{q};
        end
    else
        
    end
end
% Build Matrices
% --------------
% Volume Matrices
M = [];
K  = zeros(ntot);
G  = cell(dim,1);
for d=1:dim
    G{d} = zeros(ntot);
end
IV = [];
% Mass Matrix
if v_flags(1)
    M = bvals_v'*(bvals_v.*(qw_V*zt));
end
if lump_bool, M = diag(sum(M)); end
% Stiffness Matrix
if v_flags(2)
    for q=1:nq_V
        bg = bgrads_v(:,:,q);
        K = K + qw_V(q) * (bg*bg');
    end
end
% Gradient Matrix
if v_flags(3)
    for q=1:nq_V
        bt = bvals_v(q,:);
        bg = bgrads_v(:,:,q);
        for d=1:dim
            G{d} = G{d} + qw_V(q) * (bg(:,d)*bt)';
        end
    end
end
% Surface Matrices
MM = cell(nf,1);
G2 = cell(nf,1);
F  = cell(nf, 1);
for f=1:nf
    ffdd = f_dofs(f,:);
    bv = bvals_s{f};
    tmf = bv'*(bv.*(qw_S{f}*zt));
    MM{f} = tmf(ffdd,ffdd);
    if s_flags(2)
        G2{f} = cell(dim,1);
        for d=1:dim
            bg = squeeze(bgrads_s{f}(:,d,:));
            G2{f}{d} = bg*(bv.*(qw_S{f}*zt));
        end
    end
end

% Assign Output Arguments
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G, IV};
% Surface Matrices
varargout{2} = {MM, G2, F};
% Quadrature Structure - Volume
varargout{3} = {qx_V, qw_V, bvals_v, bgrads_v};
% Quadrature Structure - Surface
varargout{4} = {qx_S, qw_S, bvals_s, bgrads_s};
% ------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxilliary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_num_dofs(dim, ord)
out = (ord+1)^dim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_dofs(dim, ord)
if dim == 2
    if ord == 1
        out = [1,2;2,3;3,4;4,1];
    elseif ord == 2
        out = [1,2,5;2,3,6;3,4,7;4,1,8];
    elseif ord == 3
        out = [1,2,5,6;2,3,7,8;3,4,9,10;4,1,11,12];
    end
elseif dim == 3
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_ref_quadrature(dim, ord)
[tqx, tqw] = get_legendre_gauss_quad(ord); nt = length(tqw);
if dim == 1
    qx = tqx; qw = tqw;
elseif dim == 2
    qw = tqw*tqw'; qw = qw(:); qw = qw / sum(qw);
    qxx = zeros(nt, nt); qxy = zeros(nt, nt);
    for i=1:nt
        qxx(i,:) = tqx(i);
        qxy(:,i) = tqx(i);
    end
    qx = [qxx(:), qxy(:)];
elseif dim == 3
    qw = zeros(nt, nt, nt); qqw = tqw*tqw';
    qxx = zeros(nt, nt, nt); qxy = zeros(nt, nt, nt); qxz = zeros(nt, nt, nt); 
    for i=1:nt
        qxx(i,:,:) = tqx(i);
        qxy(:,i,:) = tqx(i);
        qxz(:,:,i) = tqx(i);
        qw(:,:,i)  = qqw*tqw(i);
    end
    qx = [qxx(:), qxy(:), qxz(:)];
    qw = qw(:); qw = qw / sum(qw);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trqx_S, trqw_S] = get_surface_ref_quadrature(dim, rqx_S, rqw_S, ref_x)
nrq = length(rqw_S); onrq = ones(nrq,1);
if dim == 2
    trqx_S = cell(4,1); trqw_S = cell(4,1);
    % Face 1
    for i=1:4
        ii = [i,mod(i,4)+1];
        x0 = ref_x(ii(1),:);
        dx = diff(ref_x(ii,:));
        trqx_S{i} = onrq*x0 + rqx_S*dx;
    end
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J,invJ] = evaluate_numerical_jacobian(dim, verts, x)
n = size(x,1);
J = cell(n,1); invJ = cell(n,1);
% Evaluate Jacobian
if dim == 2
    db = eval_quad_grads_for_jac(1, x);
    for i=1:n
        JJ = [db(:,1,i)'*verts(:,1), db(:,2,i)'*verts(:,1);...
                  db(:,1,i)'*verts(:,2), db(:,2,i)'*verts(:,2)];
        J{i} = JJ;
        detJ = JJ(1,1)*JJ(2,2)-JJ(2,1)*JJ(1,2);
        invJ{i} = [JJ(2,2),-JJ(1,2);-JJ(2,1),JJ(1,1)]/detJ;
    end
elseif dim == 3
    db = eval_hex_grads_for_jac(1, x);
    for i=1:n
        JJ = [db(:,1,i)'*verts(:,1), db(:,2,i)'*verts(:,1), db(:,3,i)'*verts(:,1);...
                  db(:,1,i)'*verts(:,2), db(:,2,i)'*verts(:,2), db(:,3,i)'*verts(:,2);...
                  db(:,1,i)'*verts(:,3), db(:,2,i)'*verts(:,3), db(:,3,i)'*verts(:,3)];
        J{i} = JJ;
        detJ = JJ(1,1)*(JJ(3,3)*JJ(2,2)-JJ(3,2)*JJ(2,3)) - JJ(2,1)*(JJ(3,3)*JJ(1,2)-JJ(3,2)*JJ(1,3)) + JJ(3,1)*(JJ(2,3)*JJ(1,2)-JJ(2,2)*JJ(1,3));
        invJ{i} = [(JJ(3,3)*JJ(2,2)-JJ(3,2)*JJ(2,3)),-(JJ(3,3)*JJ(1,2)-JJ(3,2)*JJ(1,3)), (JJ(2,3)*JJ(1,2)-JJ(2,2)*JJ(1,3));...
               -(JJ(3,3)*JJ(2,1)-JJ(3,1)*JJ(2,3)), (JJ(3,3)*JJ(1,1)-JJ(3,1)*JJ(1,3)),-(JJ(2,3)*JJ(1,1)-JJ(2,1)*JJ(1,3));...
                (JJ(3,2)*JJ(2,1)-JJ(3,1)*JJ(2,2)),-(JJ(3,2)*JJ(1,1)-JJ(3,1)*JJ(1,2)), (JJ(2,2)*JJ(1,1)-JJ(2,1)*JJ(1,2))]/detJ;
    end
end
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
    out(:,:,i) = [t-1,s-1;1-t,-s;t,s;-t,1-s];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%