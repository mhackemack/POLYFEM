%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D Lagrange Basis Function Generator - Volume Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to determine the number of edges for a
%                   geometric cell.
%
%   Notes:          1D -   E = 2
%                   2D -   E = V
%                   3D -   E = V + F - 2 (Euler's Polyhedra Formula)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = number_of_edges(verts, faces)
% Get Number of Vertices
% ----------------------
[m, n] = size(verts);
if n > m, verts = verts'; end
[nv, dim] = size(verts);
% Get Number of Faces
% -------------------
if iscell(faces)
    nf = length(faces);
else
    nf = size(faces,1);
end
% Set Output
% ----------
if dim == 1
    out = 2;
elseif dim == 2
    out = nv;
elseif dim == 3
    out = nv + nf - 2;
end