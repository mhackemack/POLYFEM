%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          PWQ Main Generation Function
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB script to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the PWQ basis functions.
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
function varargout = bf_cell_func_PWQ( varargin )
% Collect Input/Output Arguments
% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
ntot = get_num_serendipity_points( dim, nverts, length(faces), order);
% Quick Error Checking
% ------------------------------------------------------------------------------
if order ~= 2, error('PWQ only defined for 2nd order.'); end
if dim ~= 2, error('PWQ only defined for 2D at this time.'); end
if nv ~= nverts, error('Something is not properly defined.'); end
% Allocate Matrix Space
% ------------------------------------------------------------------------------
znv = zeros(ntot);
M = znv;
K = znv;
G = cell(dim, 1);
for d=1:dim, G{d} = znv; end
MM = cell(nf, 1);
G2 = cell(nf, 1);
for f=1:nf
    MM{f} = zeros(length(faces{f}) + (order-1));   % 2D only...
    for d=1:dim, G2{f}{d} = znv; end
end
% Collect all Matrices and Quadratures
% ------------------------------------------------------------------------------
f_dofs = get_face_dofs(nv, faces, order);
% Calculate quadrature points and basis function values
[qx_v, qw_v] = get_general_volume_quadrature(verts, faces, q_ord, true); nqx = length(qw_v);
[bmv, gmv] = PWLD_O2_basis_functions(verts, qx_v, faces, order, nverts);
% Build surface quadrature
qx_s = cell(nf,1); qw_s = cell(nf,1); qs_ind = cell(nf,1);
qxs_list = []; qws_list = []; counter = 0;
for f=1:nf
    [qx_s{f}, qw_s{f}] = get_general_surface_quadrature(verts, faces{f}, q_ord);
    qxs_list = [qxs_list; qx_s{f}]; qws_list = [qws_list; qw_s{f}];
    nx = length(qw_s{f}); tind = 1:nx; qs_ind{f} = counter + tind;
    counter = counter + nx;
end
[tbms, tgms] = PWLD_O2_basis_functions(verts, qxs_list, faces, order, nverts);
bms = cell(nf,1); gms = cell(nf,1);
for f=1:nf
    bms{f} = tbms(qs_ind{f},:);
    gms{f} = tgms(:,:,qs_ind{f});
end
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
% surface matrices
for f=1:nf
    nqx = length(qw_s{f}); fv = f_dofs{f};
    tbmsf = bms{f};
    for q=1:nqx
        bt = tbmsf(q,fv);
        MM{f} = MM{f} + qw_s{f}(q) * (bt'*bt);
        if s_flags(2)
            gt = gms{f};
            for d=1:dim
                G2{f}{d}(:,fv) = G2{f}{d}(:,fv) + qw_s{f}(q) * gt(:,d,q)*bt;
            end
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