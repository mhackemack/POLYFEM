%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to produce the basis function values and
%                   gradients for the linear case.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PWLD_O1_basis_functions(varargin)
% Collect Input Arguments
% ------------------------------------------------------------------------------
verts  = varargin{1};
qx     = varargin{2};
faces  = varargin{3};
order  = varargin{4};
nverts = varargin{5};
% Determine Input Characteristics and Perform Error Checking
% ------------------------------------------------------------------------------
[nv, dim] = size(verts); nf = length(faces);
nqx = size(qx, 1);
if order > 1, error('Only 1st order in this functor. Go fix your code.'); end
if nv ~= nverts, error('Number of vertices does not align. Go fix your code.'); end
% Allocate Memory Space
% ------------------------------------------------------------------------------
bout = zeros(nqx, nv); gout = zeros(nv, dim, nqx);
% Retrieve 1D Values and Gradients
% ------------------------------------------------------------------------------
if dim == 1
    x0 = verts(1);
    J = get_simplex_jacobian(verts,dim); invJ = 1/J;
    
end
% Retrieve 2D Values and Gradients
% ------------------------------------------------------------------------------
if dim == 2
    
end
% Retrieve 3D Values and Gradients
% ------------------------------------------------------------------------------
if dim == 3
    
end
% Assign Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if nargout > 1, varargout{2} = gout; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = get_simplex_jacobian(verts, dim)
J = zeros(dim);
vverts = verts';
for d=1:dim
    J(:,d) = vverts(:,d+1) - vverts(:,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_values(dim, xin)
if dim == 1
    out = [1-xin,xin];
elseif dim == 2
    x = xin(:,1); y = xin(:,2);
    out = [1-x-y, x, y];
else
    x = xin(:,1); y = xin(:,2); z = xin(:,3);
    out = [1-x-y-z, x, y, z];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_grads(dim)
if dim == 1
    out = [-1;1];
elseif dim == 2
    out = [-1,-1;1,0;0,1];
else
    out = [-1,-1,-1;1,0,0;0,1,0;0,0,1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%