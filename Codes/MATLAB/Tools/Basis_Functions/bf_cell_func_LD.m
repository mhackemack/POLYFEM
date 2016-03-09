%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          LD Main Generation Function
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the LD basis functions.
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
function varargout = bf_cell_func_LD( varargin )
% Collect Input/Output Arguments
% ------------------------------------------------------------------------------
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
% Quick Error Checking
% ------------------------------------------------------------------------------
if order ~= 1, error('LD requires 1st order FE space.'); end
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
[xmean,dx] = get_geometry_components(dim,verts);
% Allocate Matrices and Quadrature Points
% ------------------------------------------------------------------------------
M = zeros(dim+1);
K = zeros(dim+1);
G = cell(dim, 1);
for d=1:dim, G{d} = zeros(dim+1); end
MM = cell(nf, 1);
G2 = cell(nf, 1);
qx_s = cell(nf, 1);
qw_s = cell(nf, 1);
bms  = cell(nf, 1);
gms  = cell(nf, 1);
for f=1:nf
    MM{f} = zeros(dim+1);
    for d=1:dim, G2{f}{d} = zeros(dim+1); end
end
% Collect all Matrices and Quadratures
% ------------------------------------------------------------------------------
% Cell-Wise Values
[qx_v, qw_v] = get_general_volume_quadrature(verts, faces, q_ord, true); nqx = length(qw_v);
[bmv,gmv] = get_LD_basis(qx_v,xmean,dx);
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
for f=1:nf
    [qx_s{f}, qw_s{f}] = get_general_surface_quadrature(verts, faces{f}, q_ord, true);
    [bms{f}, gms{f}] = get_LD_basis(qx_s{f},xmean,dx);
    nqx = length(qw_s{f});
    for q=1:nqx
        bt = bms{f}(q,:);
        MM{f} = MM{f} + qw_s{f}(q) * (bt'*bt);
    end
    if s_flags(2)
        for q=1:nqx
            bt = bms{f}(q,:);
            for d=1:dim
                gt = gms{f}(:,d,q);
                G2{f}{d} = G2{f}{d} + qw_s{f}(q) * gt*bt;
            end
        end
    end
end
% Set Output Structures
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G};
% Surface Matrices
varargout{2} = {MM, G2};
% Quadrature Structures
varargout{3} = {qx_v, qw_v, bmv, gmv};
varargout{4} = {qx_s, qw_s, bms, gms};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = compute_volume_centroid(dim, verts, faces)
if dim == 1
    out = (verts(2) - verts(1))/2;
elseif dim == 2
    [out,~] = polygonCentroid(verts);
elseif dim == 3
    out = polyhedronCentroid(verts, faces);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = compute_face_centroid(dim, fverts)
if dim == 1
    out = fverts;
elseif dim == 2
    out = mean(verts);
elseif dim == 3
    out = polygonCentroid3d(fverts);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xmean,dx] = get_geometry_components(dim,verts)
xmean = mean(verts);
% xmean = compute_volume_centroid(dim, verts, faces);
if dim == 1
    dx = verts(2) - verts(1);
elseif dim == 2
    dx = [max(verts(:,1)) - min(verts(:,1)),...
          max(verts(:,2)) - min(verts(:,2))];
elseif dim == 3
    dx = [max(verts(:,1)) - min(verts(:,1)),...
          max(verts(:,2)) - min(verts(:,2)),...
          max(verts(:,3)) - min(verts(:,3))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bout, gout] = get_LD_basis(x,xmean,dx)
[nx,dim] = size(x); onx = ones(nx,1);
ddr = xmean./dx;
bout = [onx, 2*x./(onx*dx) - 2*onx*ddr];
gout = zeros(dim+1,dim,nx);
gt = [zeros(1,dim);diag(2./dx)];
for q=1:nx
    gout(:,:,q) = gt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%