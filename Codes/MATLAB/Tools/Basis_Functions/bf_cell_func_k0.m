%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          k=0 Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
%
%   Description:    MATLAB script to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the 0th-order basis functions.
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
function varargout = bf_cell_func_k0( varargin )
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
% Quick Error Checking
% ------------------------------------------------------------------------------
if ord ~= 0, error('Only 0th-order basis functions.'); end
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
h = get_max_diamter( verts );
% ------------------------------------------------------------------------------
% Allocate Matrix Space
% ------------------------------------------------------------------------------
M = 0;
K = 0;
G = cell(dim, 1);
for d=1:dim, G{d} = 0; end
MM = cell(nf, 1);
G2 = cell(nf, 1);
for f=1:nf
    MM{f} = 0;
    for d=1:dim, G2{f}{d} = 0; end
end
% Collect all Matrices and Quadratures
% ------------------------------------------------------------------------------
% Cell-Wise Values
[qx_v, qw_v] = get_general_volume_quadrature(verts, faces, q_ord, true); nqx = length(qw_v);
% mass matrix
M = sum(qw_v);
% Face-Wise Values
% ------------------------------------------------------------------------------
[qx_s, qw_s, bms, gms] = get_surface_values(dim, verts, faces, q_ord, q_ord);
for f=1:nf
    MM{f} = sum(qw_s{f});
    if s_flags(2)
        for d=1:dim
            G2{f}{d} = 0;
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
function [qx_s, qw_s, bms, gms] = get_surface_values(dim, verts, faces, ord, q_ord)
nf = length(faces);
qx_s = cell(nf, 1);
qw_s = cell(nf, 1);
bms  = cell(nf, 1);
gms  = cell(nf, 1);
% Change based on dimension type
if dim == 1
    
elseif dim == 2
    
elseif dim == 3
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%