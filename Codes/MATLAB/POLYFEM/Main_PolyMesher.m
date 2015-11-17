%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          PolyMesh Run Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate regular polygonal meshes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc; fclose all;
dim = 2; out_dir = 'geometry_inputs/precompiled/';
L = 1; dom = @SqDomain_L1; gname = 'PolyMesh_SqDomain';
tol = 1e-5; maxit = 2e3;
nele = [256^2];
% nele = [2^2, 4^2, 8^2, 16^2, 32^2];
% nele = [64^2, 128^2, 256^2];

% Loop through number of elements list
for i=1:length(nele)
    n = nele(i);
    % Generate PolyMesh
    [N,E,S,~,~] = PolyMesher(dom, n, maxit, tol);
    % Correct Data Points
    nn = size(N,1);
    for j=1:nn
        % xmin
        if abs(N(j,1) - 0) < 1e-8
            N(j,1) = 0;
        end
        % xmax
        if abs(N(j,1) - L) < 1e-8
            N(j,1) = L;
        end
        % ymin
        if abs(N(j,2) - 0) < 1e-8
            N(j,2) = 0;
        end
        % ymax
        if abs(N(j,2) - L) < 1e-8
            N(j,2) = L;
        end
    end
    % Build General Geometry Mesh
    geometry = GeneralGeometry(dim, 'polymesh', N, E);
    % Save Mesh
    save([out_dir,sprintf('%s_L%d_n%d',gname,L,n),'.mat'], 'geometry');
    clear geometry;
end