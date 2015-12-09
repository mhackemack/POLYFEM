%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Generalized Gaussian Quadrature on Polygons
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_general_polygonal_gauss_quadrature(verts, faces, deg)
% Get General Input Information 
% ------------------------------------------------------------------------------
[nverts, dim] = size(verts);
nf = length(faces);
% Quick error and problem checking
% ------------------------------------------------------------------------------
if dim ~= 2,     error('Method is only applicable for 2D polygons.'); end
if nverts ~= nf, error('Number vertices != number faces.'); end

% ------------------------------------------------------------------------------