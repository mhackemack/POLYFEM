%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Maximum Entropy Main Generation Function
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the Lagrange basis functions.
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
function varargout = bf_cell_func_max_entropy( varargin )
% Collect Input/Output Arguments
% ------------------------------------------------------------------------------
nout = nargout;
nverts = varargin{1};
verts = varargin{2}(1:nverts,:);
faces = varargin{3}; nf = length(faces);
order = varargin{4};
lump_bool = varargin{5};
v_flags = varargin{6};
s_flags = varargin{7};
q_bool = varargin{8};
q_ord = order+2;
if nargin > 8
    if ~isempty(varargin{9}),q_ord = varargin{9};end
end
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
ntot = get_num_serendipity_points( dim, nverts, length(faces), order);
% Quick Error Checking
% ------------------------------------------------------------------------------
if order > 2 && dim == 2, error('Maximum Entropy only defined for order 1 and 2 in 2D.'); end
if order > 1 && dim == 3, error('Maximum Entropy only defined for order 1 in 3D.'); end
% Allocate Matrix Space
% ------------------------------------------------------------------------------
znv = zeros(ntot);
M = znv;
K = znv;
G = cell(dim, 1);
for d=1:dim, G{d} = znv; end
MM = cell(nf, 1);
G2 = cell(nf, 1);
F  = cell(nf, 1);
for f=1:nf
    if dim == 2
        MM{f} = zeros(length(faces{f}) + (order-1));
        F{f} = zeros(length(faces{f}) + (order-1));
    else
        MM{f} = zeros(length(faces{f}));
        F{f} = zeros(length(faces{f}));
    end
    for d=1:dim, G2{f}{d} = znv; end
end
% Collect all Matrices and Quadratures
% ------------------------------------------------------------------------------
h = get_max_diamter( verts );
f_dofs = get_face_dofs(nv, faces, order);
% Calculate quadrature points and basis function values
[qx_v, qw_v] = get_general_volume_quadrature(verts, faces, q_ord, true); nqx = length(qw_v);
[bmv, gmv] = max_entropy_basis_functions(verts, qx_v, faces, order, nverts);
% if order == 1
%     [bmv, gmv] = max_entropy_O1_basis_functions(verts, qx_v, faces, order, nverts);
% elseif order == 2
%     [bmv, gmv] = max_entropy_O2_basis_functions(verts, qx_v, faces, order, nverts);
% end
% Calculate surface basis functions and gradients
[qx_s, qw_s, bms, gms] = get_ME_surface_values(dim, verts, faces, order, q_ord, h, s_flags(2));
% mass matrix
for q=1:nqx
    bt = bmv(q,:);
    M = M + qw_v(q) * (bt'*bt);
end
if lump_bool, M = diag(sum(M)); end
% stiffness matrix
if v_flags(2)
    for q=1:nqx
        bg = gmv(:,:,q);
        K = K + qw_v(q) * (bg*bg');
    end
end
% gradient matrix
if v_flags(3)
    for q=1:nqx
        bt = bmv(q,:);
        bg = gmv(:,:,q);
        for d=1:dim
            G{d} = G{d} + qw_v(q) * (bg(:,d)*bt)';
        end
    end
end
% surface matrices
for f=1:nf
    nqx = length(qw_s{f});
    fv = f_dofs{f};
    for q=1:nqx
        bt = bms{f}(q,:);
        MM{f} = MM{f} + qw_s{f}(q) * (bt'*bt);
        if s_flags(2)
            gt = gms{f};
            for d=1:dim
                G2{f}{d}(:,fv) = G2{f}{d}(:,fv) + qw_s{f}(q) * gt(:,d,q)*bt;
            end
        end
    end
    % Modify boundary values
    tbms = bms{f};
    bms{f} = zeros(nqx, ntot);
    bms{f}(:,fv) = tbms;
end
% Process Output Structures
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G};
% Surface Matrices
varargout{2} = {MM, G2, F};
% Quadrature Structures
varargout{3} = {qx_v, qw_v, bmv, gmv};
varargout{4} = {qx_s, qw_s, bms, gms};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxiallary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_max = get_max_diamter( verts )
nv = size(verts,1);
out_max = 0;
for i=1:nv
    vi = verts(i,:);
    for j=1:nv
        if i==j, continue; end
        h = norm(verts(j,:) - vi);
        if h > out_max, out_max = h; end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_dofs(nv, faces, ord)
if ord == 1
    out = faces;
else % only 2D allowed here
    nf = length(faces);
    out = cell(nf,1);
    for f=1:nf
        out{f} = [faces{f},nv+f];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx_s, qw_s, bms, gms] = get_ME_surface_values(dim, verts, faces, ord, q_ord, h, sgrad_bool)
nf = length(faces);
qx_s = cell(nf, 1);
qw_s = cell(nf, 1);
bms  = cell(nf, 1);
gms  = cell(nf, 1);
if dim == 1
    h = verts(2) - verts(1);
    qx_s{1} = verts(1); qx_s{2} = verts(2);
    qw_s{1} = 1; qw_s{2} = 1;
    bms{1} = 1; bms{2} = 1;
    if sgrad_bool, gms{1} = [-1/h;1/h]; gms{2} = [1/h;-1/h]; end
elseif dim == 2
    [tqx, tqw] = get_legendre_gauss_quad(q_ord); ntqx = length(tqw);
    ttqx1 = []; ttqx2 = [];
    fones = ones(ntqx,1);
    % Loop through faces and build some information
    for f=1:nf
        fv = faces{f};
        v = verts(fv,:);
        dx = v(2,:) - v(1,:);
        len = norm(diff(v));
        n = [dx(2), -dx(1)]/len;
        qw_s{f} = tqw*len;
        qx_s{f} = fones*v(1,:) + tqx*dx;
        ttqx1 = [ttqx1;qx_s{f} - fones*n*h/1e2];
%         ttqx2 = [ttqx2;qx_s{f} - fones*n*h/2e2];
        if ord == 1
            bms{f} = [1-tqx, tqx];
        elseif ord == 2
             bms{f} = [(1-tqx).^2, tqx.^2, 2*tqx.*(1-tqx)];
        end
    end
    % Get Gradient Estimates
    if sgrad_bool
        if ord == 1
            [~,tg] = max_entropy_O1_basis_functions(verts, ttqx1, faces, ord, size(verts,1));
        elseif ord == 2
            [~,tg] = max_entropy_O2_basis_functions(verts, ttqx1, faces, ord, size(verts,1));
        end
        % Rebuild Surface Gradients
        for f=1:nf
            iif = ntqx*(f-1)+1:ntqx*f;
            gms{f} = tg(:,:,iif);
        end
    end
elseif dim == 3
    for f=1:nf
        fv = faces{f};
        v = verts(fv,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%