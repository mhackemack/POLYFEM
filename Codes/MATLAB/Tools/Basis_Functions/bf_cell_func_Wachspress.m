%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Wachspress Main Generation Function
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
%                   5) Volumetric Matrix Flags
%                   6) Surface Matrix Flags
%                   7) Quadrature boolean
%                   8) Quadrature Order (Optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = bf_cell_func_Wachspress( varargin )
% Collect Input/Output Arguments
% ------------------------------
nout = nargout;
nverts = varargin{1};
verts = varargin{2}(1:nverts,:);
faces = varargin{3}; nf = length(faces);
order = varargin{4};
v_flags = varargin{5};
s_flags = varargin{6};
q_bool = varargin{7};
q_ord = order+2;
if nargin > 7
    if ~isempty(varargin{8}),q_ord = varargin{8};end
end
% Prepare Vertices and Dimensional Space
% --------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
% Quick Error Checking
% --------------------
if order > 1, error('Wachpress only defined for order 1.'); end
% Compute and exit immediately if 1D
% ------------------------------------------------------------------------------
if dim == 1
    [bf_V,bf_S,QV,QS] = bf_cell_func_1D(varargin{:});
    varargout = {bf_V, bf_S, QV, QS};
    return
end
% ------------------------------------------------------------------------------
% Allocate Matrix Space
% ---------------------
znv = zeros(nv);
J = zeros(dim);
M = znv;
K = znv;
G = cell(dim, 1);
for d=1:dim, G{d} = znv; end
MM = cell(nf, 1);
G2 = cell(nf, 1);
for f=1:nf
    MM{f} = zeros(length(faces{f}));
    for d=1:dim, G2{f}{d} = znv; end
end
% Collect all Matrices and Quadratures
% ------------------------------------------------------------------------------
% Cell-Wise Values
[qx_v, qw_v] = get_general_volume_quadrature(verts, faces, 2*order+1, true); nqx = length(qw_v);
[bmv, gmv] = wachspress_basis_functions(verts, qx_v, faces, order, nverts);
% mass matrix
for q=1:nqx
    bt = bmv(q,:);
    M = M + qw_v(q) * (bt'*bt);
end
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
% Face-Wise Values
[qx_s, qw_s, bms, gms] = get_surface_values(dim, verts, faces, order);
for f=1:nf
    nqx = length(qw_s{f});
    for q=1:nqx
        bt = bms{f}(q,:);
        MM{f} = MM{f} + qw_s{f}(q) * (bt'*bt);
        if s_flags(2)
            
        end
    end
end
% Process Output Structures
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G};
% Surface Matrices
varargout{2} = {MM, G2};
% Quadrature Structures
varargout{3} = {qx_v, qw_v, bmv, gmv};
varargout{4} = {qx_s, qw_s, bms, gms};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxiallary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx_s, qw_s, bms, gms] = get_surface_values(dim, verts, faces, ord)
nf = length(faces);
qx_s = cell(nf, 1);
qw_s = cell(nf, 1);
bms  = cell(nf, 1);
gms  = cell(nf, 1);
if dim == 2
    [tqx, tqw] = lgwt(ord+2,0,1); ntqx = length(tqw);
    fones = ones(ntqx,1);
    for f=1:nf
        fv = faces{f};
        v = verts(fv,:);
        dv = v(2,:) - v(1,:);
        len = norm(diff(v));
        qw_s{f} = tqw*len;
        qx_s{f} = fones*v(1,:) + tqx*dv;
        if ord == 1
            bms{f} = [1-tqx, tqx];
        elseif ord == 2
             bms{f} = [(1-tqx).^2, tqx.^2, 2*tqx.*(1-tqx)];
        end
    end
elseif dim == 3
    for f=1:nf
        fv = faces{f};
        v = verts(fv,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%