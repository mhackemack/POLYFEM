%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Harmonic Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the basis function and gradient
%                   values using the Maximum Entropy Methodology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = harmonic_basis_functions( varargin )
nout = nargout;
grad_bool = false;
% Collect Input Arguments
% -----------------------
verts = varargin{1};
qx = varargin{2};
faces = varargin{3};
% Prepare Vertices and Dimensional Space
% --------------------------------------
[nverts, dim] = size(verts);
nqx = size(qx, 1);
% Allocate Matrix Memory
% ----------------------
if nout > 1, grad_bool = true; end
ntot = get_num_serendipity_points( dim, nverts, length(faces), 1);
bout = zeros(nqx, ntot);
if grad_bool, gout = zeros(ntot, dim, nqx); end
% Get Problem Preliminaries
% -------------------------
[rv, sv] = get_vertex_differences(verts, qx);
% Build Basis Function Sets
% -------------------------

if grad_bool, gout = get_basis_gradients(bout, R); end
% Set Output Arguments
% --------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rv, sv] = get_vertex_differences( verts, qx )
[nv, dim] = size(verts);
nqx = size(qx, 1);
rv = zeros(nv, dim, nqx);
sv = zeros(nv, nqx);
zn = ones(nv, 1);
zones = ones(dim,1);
for q=1:nqx
    rv(:,:,q) = verts - zn*qx(q,:);
    rrv = rv(:,:,q).*rv(:,:,q);
    sv(:,q) = sqrt(rrv*zones);
end
sv = sv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%