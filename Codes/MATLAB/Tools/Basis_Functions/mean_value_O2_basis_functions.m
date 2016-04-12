%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Quadratic Mean Value Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the basis function and gradient
%                   values using the quadratic mean value methodology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = mean_value_O2_basis_functions(varargin)
% Collect Input Arguments
% ------------------------------------------------------------------------------
verts  = varargin{1};
qx     = varargin{2};
faces  = varargin{3};
order  = varargin{4};
nverts = varargin{5};
% Determine Input Characteristics
% ------------------------------------------------------------------------------
[nv, dim] = size(verts);
if dim ~= 2, error('2nd order MV is only for 2D problems at this time.'); end
if order ~= 2, error('Only 2nd order in this functor. Go fix your code.'); end
if nv ~= nverts, error('Number of vertices does not align. Go fix your code.'); end
% Call 
% ------------------------------------------------------------------------------
nout = nargout;
if nout == 1
    bout = barycentric_serendipity_Rev2(verts, qx, faces, @mean_value_O1_basis_functions);
elseif nout == 2
    [bout, gout] = barycentric_serendipity_Rev2(verts, qx, faces, @mean_value_O1_basis_functions);
end
% Assign Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if nout > 1, varargout{2} = gout; end