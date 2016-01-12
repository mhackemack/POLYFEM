%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D Basis Function Plotter Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Project Space
% ------------------------------------------------------------------------------
if exist('pbool', 'var')
    clearvars -except pbool
else
    clear; pbool = false;
end
clc; close all; format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Populate global space
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Office');
% ------------------------------------------------------------------------------
% Begin user input section
% ------------------------------------------------------------------------------
% Input/Output informatin
% -----------------------
out_dir = 'outputs/BF_Plots/';
% Basis function information
% --------------------------
p = 1;
% BF1
BF(p).Name = 'PWLD';
BF(p).Degree = 1;
BF(p).BasisFunction = @PWLD_basis_functions;
p = p + 1;
% BF2
BF(p).Name = 'PWLD';
BF(p).Degree = 2;
BF(p).BasisFunction = @PWLD_basis_functions;
p = p + 1;
% BF3
BF(p).Name = 'WACHSPRESS';
BF(p).Degree = 1;
BF(p).BasisFunction = @wachspress_basis_functions;
p = p + 1;
% BF4
BF(p).Name = 'WACHSPRESS';
BF(p).Degree = 2;
BF(p).BasisFunction = @wachspress_basis_functions;
p = p + 1;
% BF5
BF(p).Name = 'MV';
BF(p).Degree = 1;
BF(p).BasisFunction = @mean_value_basis_functions;
p = p + 1;
% BF6
BF(p).Name = 'MV';
BF(p).Degree = 2;
BF(p).BasisFunction = @mean_value_basis_functions;
p = p + 1;
% BF7
BF(p).Name = 'MAXENT';
BF(p).Degree = 1;
BF(p).BasisFunction = @max_entropy_basis_functions;
p = p + 1;
% BF8
BF(p).Name = 'MAXENT';
BF(p).Degree = 2;
BF(p).BasisFunction = @max_entropy_basis_functions;
% % BF9
% BF(p).Name = 'LAGRANGE';
% BF(p).Degree = 1;
% BF(p).BasisFunction = get_lagrange_function('vals',2,2);
% p = p + 1;
% % BF10
% BF(p).Name = 'LAGRANGE';
% BF(p).Degree = 2;
% BF(p).BasisFunction = get_lagrange_function('vals',2,2);
% p = p + 1;
% Geometric cell information
% --------------------------
% Square quadrature
[xs,ys] = meshgrid(0,1,101);
xxs = xs(:); yys = ys(:);
qxs = [xxs,yys];
p = 1;
% Geometry1
Geometry(p).Name = 'square';
Geometry(p).Vertices = [0,0;1,0;1,1;0,1];
Geometry(p).Faces = {[1,2],[2,3],[3,4],[4,1]};
Geometry(p).Quad = qxs;
p = p + 1;
% Geometry2
Geometry(p).Name = 'deg_square';
Geometry(p).Vertices = [0,0;1,0;1,1;.5,1;0,1];
Geometry(p).Faces = {[1,2],[2,3],[3,4],[4,5],[5,1]};
Geometry(p).Quad = qxs;
p = p + 1;
% Geometry3
Geometry(p).Name = 'L-domain';
Geometry(p).Vertices = [0,0;1,0;1,0.5;0.5,0.5;0.5,1;0,1];
Geometry(p).Faces = {[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]};
Geometry(p).Quad = qxs;
p = p + 1;
% ------------------------------------------------------------------------------
% End user input section
% ------------------------------------------------------------------------------
% Perform all plotting procedures
% ------------------------------------------------------------------------------
% Gather some basic information
ng = length(Geometry);
nbf = length(BF);
% Loop through geometries
for g=1:ng
    geo = Geometry(g);
    % Determine if quadrature points are inside the domain
    inp = inpolygon(geo.Quad(:,1), geo.Quad(:,2), geo.Vertices(:,1), geo.Vertices(:,2));
    % Loop through basis functions
    for b = 1:nbf
        % Generate basis function values
        
        % Generate final output name
        out_name = [geo.Name,'_',BF(b).name,num2str(BF(b).Degree)];
        
        % Generate and save all output
        
    end
end