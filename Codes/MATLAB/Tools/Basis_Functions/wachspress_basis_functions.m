%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Wachpress Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the basis function and gradient
%                   values using the Wachpress Methodology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = wachspress_basis_functions( varargin )
% Collect Input Arguments
% ------------------------------------------------------------------------------
verts  = varargin{1};
qx     = varargin{2};
faces  = varargin{3};
order  = varargin{4};
nverts = varargin{5};
grad_bool = false;
if nargout > 1, grad_bool = true; end
% Get Basis Function Values and Gradients from Appropriate Functor Calls
% ------------------------------------------------------------------------------
if order == 1 
    if ~grad_bool
        bout = wachspress_O1_basis_functions(verts, qx, faces, order, nverts);
    else
        [bout, gout] = wachspress_O1_basis_functions(verts, qx, faces, order, nverts);
    end
elseif order == 2
    if ~grad_bool
        bout = wachspress_O2_basis_functions(verts, qx, faces, order, nverts);
    else
        [bout, gout] = wachspress_O2_basis_functions(verts, qx, faces, order, nverts);
    end
end
% Set Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end