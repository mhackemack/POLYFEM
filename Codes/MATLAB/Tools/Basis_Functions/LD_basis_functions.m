%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          LD Basis Function Generator
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
function varargout = LD_basis_functions(varargin)
% Collect Input Arguments
% ------------------------------------------------------------------------------
verts  = varargin{1};
qx     = varargin{2};
faces  = varargin{3};
order  = varargin{4};
nverts = varargin{5};
grad_bool = false;
if nargout > 1, grad_bool = true; end

% Allocate Memory Space
% ------------------------------------------------------------------------------
[nqx,dim] = size(qx);
[xmean,dx] = get_geometry_components(dim,verts);
% Set Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end


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