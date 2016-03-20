%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          PolyMesh (Uniform) Run Script
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
if exist('pbool', 'var')
    clearvars -except pbool
else
    clear; pbool = false;
end
clc; close all; format long e; clear persistent;
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
dim = 2;
out_dir = 'geometry_inputs/precompiled/';
gname = 'PolyMesh_SqDomain_Uniform';
L = 1; dom = @SqDomain_L1;
nxy = [2,4,8,16];
tol = 1e-6; maxit = 1e4;
% nxy = [2,4,8,16,32,64,128,256];
% Loop through number of elements list
for i=1:length(nxy)
    n = nxy(i);
    dx = L/n; dy = L/n;
    [x,y] = meshgrid(dx/2:dx:L, dy/2:dy:L);
    P = [x(:),y(:)];
    [N,E,S,~,~] = PolyMesher(dom, n^2, maxit, tol, P);
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
    save([out_dir,sprintf('%s_L%d_n%d',gname,L,n^2),'.mat'], 'geometry');
    clear geometry;
end