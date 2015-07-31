%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Orthogonal Cartesian Lagrange Main Generation Function
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
%                   5) Volumetric Matrix Flags
%                   6) Surface Matrix Flags
%                   7) Quadrature boolean
%                   8) Quadrature Order (Optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = bf_cell_func_Lagrange_Cartesian_Orthogonal( varargin )
% Collect Input/Output Arguments
% ------------------------------
nout = nargout;
nverts = varargin{1};
verts = varargin{2}(1:nverts,:);
faces = varargin{3}; nf = length(faces);
ord = varargin{4};
v_flags = varargin{5};
s_flags = varargin{6};
q_bool = varargin{7};
q_ord = ord+2;
if nargin > 7
    if ~isempty(varargin{8}),q_ord = varargin{8};end
end
% Prepare Vertices and Dimensional Space
% --------------------------------------
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
ntot = get_num_dofs(dim, ord); zt = ones(1,ntot);
[g_to_r, r_to_g, f_g_to_r, f_r_to_g] = get_proper_indexing(dim, verts, faces, ord, ntot);
% Retrieve Reference Function Handles
% -----------------------------------
ref_val_func = get_lagrange_function('vals', 2, dim);
ref_grad_func = get_lagrange_function('grads', 2, dim);
% Generate Reference Quadrature Spaces
% ------------------------------------
if dim == 2
    ref_x = [0,0;1,0;1,1;0,1];
elseif dim == 3
    ref_x = [0,0,0;1,0,0;1,1,0;0,1,0;0,0,1;1,0,1;1,1,1;0,1,1];
end
[rqx_V, rqw_V] = get_ref_quadrature(dim, q_ord);   nq_V = length(rqw_V);
[rqx_S, rqw_S] = get_ref_quadrature(dim-1, q_ord); nq_S = length(rqw_S);
[trqx_S, trqw_S] = get_surface_ref_quadrature(dim, rqx_S, rqw_S, ref_x);
zs = ones(nq_S, 1);
% Compute Jacobian - constant for orthogonal cell
% -----------------------------------------------
x0 = verts(r_to_g(1),:);
J = zeros(dim); invJ = zeros(dim);
dx = max(verts(:,1)) - min(verts(:,1));
dy = max(verts(:,2)) - min(verts(:,2));
svol = dx*dy; detJ = svol;
J(1,1) = dx; J(2,2) = dy;
invJ(1,1) = 1/dx; invJ(2,2) = 1/dy;
% Compute Real Space Quadrature and Basis Values/Gradients
% --------------------------------------------------------
% Volume
qx_V = zeros(nq_V, dim); qw_V = rqw_V*svol;
bvals_v = ref_val_func(ord, rqx_V);
bgrads_v = zeros(ntot,dim,nq_V);
for q=1:nq_V
    qx_V(q,:) = x0 + (J*rqx_V(q,:)')';
    tg = ref_grad_func(ord, rqx_V(q,:));
    bgrads_v(:,:,q) = tg*invJ;
end
% Surface
f_dofs = get_face_dofs(dim, ord);
bvals_s = cell(nf, 1); bgrads_s = cell(nf, 1);
qx_S    = cell(nf, 1); qw_S     = cell(nf, 1);
for f=1:nf
    fv = faces{f}; fverts = verts(fv,:);
    f0 = fverts(1,:);
    bvals_s{f} = ref_val_func(ord,trqx_S{f});
    bgrads_s{f} = zeros(ntot, dim, nq_S);
    for q=1:nq_S
        bgrads_s{f}(:,:,q) = ref_grad_func(ord, trqx_S{f}(q,:))*invJ;
    end
    if dim == 2
        ddx = diff(fverts);
        len = norm(ddx);
        qw_S{f} = rqw_S*len;
        qx_S{f} = zs*f0 + rqx_S*ddx;
    else
        
    end
end
% Build Matrices
% --------------
% Volume Matrices
[M,K,G] = get_volume_matrices(dim, ord, dx, dy);
M = M(r_to_g,r_to_g);
K = K(r_to_g,r_to_g);
for d=1:dim
    G{d} = G{d}(r_to_g,r_to_g);
end
% Surface Matrices
MM = cell(nf,1);
G2 = cell(nf,1);
for f=1:nf
    bv = bvals_s{f};
    tmf = bv'*(bv.*(qw_S{f}*zt));
    MM{f} = tmf(f_r_to_g{f},f_r_to_g{f});
end
% Assign Output Arguments
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G};
% Surface Matrices
varargout{2} = {MM, G2};
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
    out = [1,2,3,4;5,8,7,6;1,5,6,2;3,7,8,4;2,6,7,3;4,8,5,1];
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
function [g_to_r, r_to_g, f_g_to_r, f_r_to_g] = get_proper_indexing(dim, verts, faces, ord, ndof)
xmin = min(verts(:,1)); xmax = max(verts(:,1));
ymin = min(verts(:,2)); ymax = max(verts(:,2));
% y
%/\
% |  4 ----- 3
% |  |       |
% |  |       |
% |  1 ----- 2
% |
% ------------> x
if dim == 2
    if size(verts,1) ~= 4, error('4 Vertices in 2D.'); end
    g_to_r = zeros(ndof,1); r_to_g = zeros(ndof,1);
    f_g_to_r = cell(4,1); f_r_to_g = cell(4,1);
    for i=1:4
        fv = faces{i};
        f_ord = 1:ord-1; nf_ord = length(f_ord);
        if     abs(xmin - verts(i,1)) < 1e-13 && abs(ymin - verts(i,2)) < 1e-13
            g_to_r(i) = 1;
            r_to_g(1) = i;
            f_g_to_r{i} = [1,2];
            f_r_to_g{1} = fv;
        elseif abs(xmax - verts(i,1)) < 1e-13 && abs(ymin - verts(i,2)) < 1e-13
            g_to_r(i) = 2;
            r_to_g(2) = i;
            f_g_to_r{i} = [2,3];
            f_r_to_g{2} = fv;
        elseif abs(xmax - verts(i,1)) < 1e-13 && abs(ymax - verts(i,2)) < 1e-13
            g_to_r(i) = 3;
            r_to_g(3) = i;
            f_g_to_r{i} = [3,4];
            f_r_to_g{3} = fv;
        elseif abs(xmin - verts(i,1)) < 1e-13 && abs(ymax - verts(i,2)) < 1e-13
            g_to_r(i) = 4;
            r_to_g(4) = i;
            f_g_to_r{i} = [4,1];
            f_r_to_g{4} = fv;
        end
    end
    % Interior Nodes
    if ord == 2
        g_to_r(9) = 9;
        r_to_g(9) = 9;
    elseif ord == 3
        rtg = r_to_g(1);
    end
else
    if size(verts,1) ~= 8, error('8 Vertices in 3D.'); end
    g_to_r = zeros(ndof,1); r_to_g = zeros(ndof,1);
    zmin = min(verts(:,3)); zmax = max(verts(:,3));
    for i=1:8
        if     abs(xmin - verts(i,1)) < 1e-13 && abs(ymin - verts(i,2)) < 1e-13 && abs(zmin - verts(i,3)) < 1e-13
            g_to_r(i) = 1;
            r_to_g(1) = i;
        elseif abs(xmax - verts(i,1)) < 1e-13 && abs(ymin - verts(i,2)) < 1e-13 && abs(zmin - verts(i,3)) < 1e-13
            g_to_r(i) = 2;
            r_to_g(2) = i;
        elseif abs(xmax - verts(i,1)) < 1e-13 && abs(ymax - verts(i,2)) < 1e-13 && abs(zmin - verts(i,3)) < 1e-13
            g_to_r(i) = 3;
            r_to_g(3) = i;
        elseif abs(xmin - verts(i,1)) < 1e-13 && abs(ymax - verts(i,2)) < 1e-13 && abs(zmin - verts(i,3)) < 1e-13
            g_to_r(i) = 4;
            r_to_g(4) = i;
        elseif abs(xmin - verts(i,1)) < 1e-13 && abs(ymin - verts(i,2)) < 1e-13 && abs(zmax - verts(i,3)) < 1e-13
            g_to_r(i) = 5;
            r_to_g(5) = i;
        elseif abs(xmax - verts(i,1)) < 1e-13 && abs(ymin - verts(i,2)) < 1e-13 && abs(zmax - verts(i,3)) < 1e-13
            g_to_r(i) = 6;
            r_to_g(6) = i;
        elseif abs(xmax - verts(i,1)) < 1e-13 && abs(ymax - verts(i,2)) < 1e-13 && abs(zmax - verts(i,3)) < 1e-13
            g_to_r(i) = 7;
            r_to_g(7) = i;
        elseif abs(xmin - verts(i,1)) < 1e-13 && abs(ymax - verts(i,2)) < 1e-13 && abs(zmax - verts(i,3)) < 1e-13
            g_to_r(i) = 8;
            r_to_g(8) = i;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,K,G] = get_volume_matrices(dim, ord, dx, dy)
if dim == 2
    if ord==1
        M = dx*dy*[ 4, 2, 1, 2;...
                    2, 4, 2, 1;...
                    1, 2, 4, 2;...
                    2, 1, 2, 4]/36;
        K = [  (dx^2 + dy^2)/(3*dx*dy),    dx/(6*dy) - dy/(3*dx), -(dx^2 + dy^2)/(6*dx*dy),    dy/(6*dx) - dx/(3*dy);...
                 dx/(6*dy) - dy/(3*dx),  (dx^2 + dy^2)/(3*dx*dy),    dy/(6*dx) - dx/(3*dy), -(dx^2 + dy^2)/(6*dx*dy);...
              -(dx^2 + dy^2)/(6*dx*dy),    dy/(6*dx) - dx/(3*dy),  (dx^2 + dy^2)/(3*dx*dy),    dx/(6*dy) - dy/(3*dx);...
                 dy/(6*dx) - dx/(3*dy), -(dx^2 + dy^2)/(6*dx*dy),    dx/(6*dy) - dy/(3*dx),  (dx^2 + dy^2)/(3*dx*dy)];
        G{1} = dy/12 * [ -2, 2, 1, -1;...
                         -2, 2, 1, -1;...
                         -1, 1, 2, -2;...
                         -1, 1, 2, -2];
        G{2} = dx/12 * [ -2, -1, 1, 2;...
                         -1, -2, 2, 1;...
                         -1, -2, 2, 1;...
                         -2, -1, 1, 2];
    elseif ord == 2
        M = dx*dy/900 * [ 16, -4,  1, -4,   8,  -2,  -2,   8,   4;...
                          -4, 16, -4,  1,   8,   8,  -2,  -2,   4;...
                           1, -4, 16, -4,  -2,   8,   8,  -2,   4;...
                          -4,  1, -4, 16,  -2,  -2,   8,   8,   4;...
                           8,  8, -2, -2,  64,   4, -16,   4,  32;...
                          -2,  8,  8, -2,   4,  64,   4, -16,  32;...
                          -2, -2,  8,  8, -16,   4,  64,   4,  32;...
                           8, -2, -2,  8,   4, -16,   4,  64,  32;...
                           4,  4,  4,  4,  32,  32,  32,  32, 256];
        K = 1/90 * [  (28*(dx^2 + dy^2))/(dx*dy),       (4*dy)/dx - (7*dx)/dy,      -(dx^2 + dy^2)/(dx*dy),       (4*dx)/dy - (7*dy)/dx,       (14*dx)/dy - (32*dy)/dx,  (2*(4*dx^2 + dy^2))/(dx*dy),   (2*(dx^2 + 4*dy^2))/(dx*dy),      (14*dy)/dx - (32*dx)/dy,   -(16*(dx^2 + dy^2))/(dx*dy);...
                           (4*dy)/dx - (7*dx)/dy,  (28*(dx^2 + dy^2))/(dx*dy),       (4*dx)/dy - (7*dy)/dx,      -(dx^2 + dy^2)/(dx*dy),       (14*dx)/dy - (32*dy)/dx,      (14*dy)/dx - (32*dx)/dy,   (2*(dx^2 + 4*dy^2))/(dx*dy),  (2*(4*dx^2 + dy^2))/(dx*dy),   -(16*(dx^2 + dy^2))/(dx*dy);...
                          -(dx^2 + dy^2)/(dx*dy),       (4*dx)/dy - (7*dy)/dx,  (28*(dx^2 + dy^2))/(dx*dy),       (4*dy)/dx - (7*dx)/dy,   (2*(dx^2 + 4*dy^2))/(dx*dy),      (14*dy)/dx - (32*dx)/dy,       (14*dx)/dy - (32*dy)/dx,  (2*(4*dx^2 + dy^2))/(dx*dy),   -(16*(dx^2 + dy^2))/(dx*dy);...
                           (4*dx)/dy - (7*dy)/dx,      -(dx^2 + dy^2)/(dx*dy),       (4*dy)/dx - (7*dx)/dy,  (28*(dx^2 + dy^2))/(dx*dy),   (2*(dx^2 + 4*dy^2))/(dx*dy),  (2*(4*dx^2 + dy^2))/(dx*dy),       (14*dx)/dy - (32*dy)/dx,      (14*dy)/dx - (32*dx)/dy,   -(16*(dx^2 + dy^2))/(dx*dy);...
                         (14*dx)/dy - (32*dy)/dx,     (14*dx)/dy - (32*dy)/dx, (2*(dx^2 + 4*dy^2))/(dx*dy), (2*(dx^2 + 4*dy^2))/(dx*dy),      (112*dx)/dy + (64*dy)/dx,  -(16*(dx^2 + dy^2))/(dx*dy),    (16*(dx^2 - dy^2))/(dx*dy),  -(16*(dx^2 + dy^2))/(dx*dy), -(32*(4*dx^2 - dy^2))/(dx*dy);...
                     (2*(4*dx^2 + dy^2))/(dx*dy),     (14*dy)/dx - (32*dx)/dy,     (14*dy)/dx - (32*dx)/dy, (2*(4*dx^2 + dy^2))/(dx*dy),   -(16*(dx^2 + dy^2))/(dx*dy),     (64*dx)/dy + (112*dy)/dx,   -(16*(dx^2 + dy^2))/(dx*dy),  -(16*(dx^2 - dy^2))/(dx*dy),  (32*(dx^2 - 4*dy^2))/(dx*dy);...
                     (2*(dx^2 + 4*dy^2))/(dx*dy), (2*(dx^2 + 4*dy^2))/(dx*dy),     (14*dx)/dy - (32*dy)/dx,     (14*dx)/dy - (32*dy)/dx,    (16*(dx^2 - dy^2))/(dx*dy),  -(16*(dx^2 + dy^2))/(dx*dy),      (112*dx)/dy + (64*dy)/dx,  -(16*(dx^2 + dy^2))/(dx*dy), -(32*(4*dx^2 - dy^2))/(dx*dy);...
                         (14*dy)/dx - (32*dx)/dy, (2*(4*dx^2 + dy^2))/(dx*dy), (2*(4*dx^2 + dy^2))/(dx*dy),     (14*dy)/dx - (32*dx)/dy,   -(16*(dx^2 + dy^2))/(dx*dy),  -(16*(dx^2 - dy^2))/(dx*dy),   -(16*(dx^2 + dy^2))/(dx*dy),     (64*dx)/dy + (112*dy)/dx,  (32*(dx^2 - 4*dy^2))/(dx*dy);...
                     -(16*(dx^2 + dy^2))/(dx*dy), -(16*(dx^2 + dy^2))/(dx*dy), -(16*(dx^2 + dy^2))/(dx*dy), -(16*(dx^2 + dy^2))/(dx*dy), -(32*(4*dx^2 - dy^2))/(dx*dy), (32*(dx^2 - 4*dy^2))/(dx*dy), -(32*(4*dx^2 - dy^2))/(dx*dy), (32*(dx^2 - 4*dy^2))/(dx*dy),   (256*(dx^2 + dy^2))/(dx*dy)];
        G{1} = dy/180 * [-12, -4,  1,   3,  16,  -2,  -4,  -6,   8;...
                           4, 12, -3,  -1, -16,   6,   4,   2,  -8;...
                          -1, -3, 12,   4,   4,   6, -16,   2,  -8;...
                           3,  1, -4, -12,  -4,  -2,  16,  -6,   8;...
                         -16, 16, -4,   4,   0,   8,   0,  -8,   0;...
                           2,  6,  6,   2,  -8,  48,  -8,  16, -64;...
                           4, -4, 16, -16,   0,   8,   0,  -8,   0;...
                          -6, -2, -2,  -6,   8, -16,   8, -48,  64;...
                          -8,  8,  8,  -8,   0,  64,   0, -64,   0];
        
        G{2} = dx/180 * [-12,   3,  1, -4,  -6,  -4,  -2,  16,   8;...
                           3, -12, -4,  1,  -6,  16,  -2,  -4,   8;...
                          -1,   4, 12, -3,   2, -16,   6,   4,  -8;...
                           4,  -1, -3, 12,   2,   4,   6, -16,  -8;...
                          -6,  -6, -2, -2, -48,   8, -16,   8,  64;...
                           4, -16, 16, -4,  -8,   0,   8,   0,   0;...
                           2,   2,  6,  6,  16,  -8,  48,  -8, -64;...
                         -16,   4, -4, 16,  -8,   0,   8,   0,   0;...
                          -8,  -8,  8,  8, -64,   0,  64,   0,   0];
    elseif ord == 3
        
    end
elseif dim == 3
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%