%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to produce the basis function values and
%                   gradients for the quadratic case.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PWLD_O2_basis_functions(varargin)
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
if dim ~= 2, error('2nd order PWLD is only for 2D problems at this time.'); end
if order ~= 2, error('Only 2nd order in this functor. Go fix your code.'); end
if nv ~= nverts, error('Number of vertices does not align. Go fix your code.'); end
% Call 
% ------------------------------------------------------------------------------
nout = nargout;
if nout == 1
    bout = barycentric_serendipity_Rev2(verts, qx, faces, @PWLD_O1_basis_functions);
elseif nout == 2
    [bout, gout] = barycentric_serendipity_Rev2(verts, qx, faces, @PWLD_O1_basis_functions);
end
% Assign Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if nout > 1, varargout{2} = gout; end