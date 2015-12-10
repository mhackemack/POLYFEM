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
function varargout = PWLD_basis_functions(varargin)
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
        bout = PWLD_O1_basis_functions(verts, qx, faces, order, nverts);
    else
        [bout, gout] = PWLD_O1_basis_functions(verts, qx, faces, order, nverts);
    end
elseif order == 2
    if ~grad_bool
        bout = PWLD_O2_basis_functions(verts, qx, faces, order, nverts);
    else
        [bout, gout] = PWLD_O2_basis_functions(verts, qx, faces, order, nverts);
    end
end
% Set Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end