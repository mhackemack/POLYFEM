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
clc; close all; format short;
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
p = p + 1;
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
ns = 301;
[xs,ys] = meshgrid(linspace(0,1,ns));
xxs = xs(:); yys = ys(:);
qxs = [xxs,yys];
p = 1;
% Geometry1
Geometry(p).Name = 'square';
Geometry(p).Vertices = [0,0;1,0;1,1;0,1];
Geometry(p).Faces = {[1,2],[2,3],[3,4],[4,1]};
Geometry(p).Quad = qxs;
Geometry(p).QuadGrid = {xs, ys};
Geometry(p).QuadPerDim = ns;
p = p + 1;
% Geometry2
Geometry(p).Name = 'deg_square';
Geometry(p).Vertices = [0,0;1,0;1,1;.5,1;0,1];
Geometry(p).Faces = {[1,2],[2,3],[3,4],[4,5],[5,1]};
Geometry(p).Quad = qxs;
Geometry(p).QuadGrid = {xs, ys};
Geometry(p).QuadPerDim = ns;
p = p + 1;
% % Geometry3
% Geometry(p).Name = 'L-domain';
% Geometry(p).Vertices = [0,0;1,0;1,0.5;0.5,0.5;0.5,1;0,1];
% Geometry(p).Faces = {[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]};
% Geometry(p).Quad = qxs;
% Geometry(p).QuadGrid = {xs, ys};
% Geometry(p).QuadPerDim = ns;
% p = p + 1;
% % Geometry4
% Geometry(p).Name = 'triangle';
% Geometry(p).Vertices = [0,0;1,0;0,1];
% Geometry(p).Faces = {[1,2],[2,3],[3,1]};
% Geometry(p).Quad = qxs;
% Geometry(p).QuadGrid = {xs, ys};
% Geometry(p).QuadPerDim = ns;
% p = p + 1;
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
    geom = Geometry(g);
    v = geom.Vertices; nv = length(v);
    f = geom.Faces;
    q = geom.Quad; qpd = geom.QuadPerDim;
    qg = geom.QuadGrid;
    % Determine if quadrature points are inside the domain
    inp = inpolygon(q(:,1), q(:,2), geom.Vertices(:,1), geom.Vertices(:,2));
    % Loop through basis functions
    for b = 1:nbf
        % Generate basis function values
        if strcmpi(BF(b).Name, 'lagrange')
            bvals = BF(b).BasisFunction(BF(b).Degree, q);
        else
            bvals = BF(b).BasisFunction(v,q,f,BF(b).Degree,nv);
        end
        bvals(~inp,:) = nan;
        % Generate final output name
        out_name = [geom.Name,'_',BF(b).Name,num2str(BF(b).Degree)];
        % Generate and save individual basis function output
        for i=1:size(bvals,2)
            bb = bvals(:,i); bbr = reshape(bb,qpd,qpd);
            % Generate contour plot
            clf; hold on; ax = gca;
            patch(v(:,1),v(:,2),[1,1,1])
            contour(qg{1}, qg{2}, bbr, 60);
            ax.XTick = (0:.1:1);
            ax.YTick = (0:.1:1);
            axis square;
            % Save contour plot
            saveas(gcf, [out_dir,out_name,'_contour_b',num2str(i)], 'fig');
            print(gcf,'-dpng',[out_dir,out_name,'_contour_b',num2str(i)]);
            print(gcf,'-depsc',[out_dir,out_name,'_contour_b',num2str(i)]);
            % Generate surf plot
            clf; hold on; ax = gca;
            surf(qg{1}, qg{2}, bbr,'linestyle','none');
%             patch(v(:,1),v(:,2),-1*ones(nv,1),[1,1,1])
            ax.XTick = (0:.1:1);
            ax.YTick = (0:.1:1);
            axis square;
            % Save surf plot
            saveas(gcf, [out_dir,out_name,'_surf_b',num2str(i)], 'fig');
            print(gcf,'-dpng',[out_dir,out_name,'_surf_b',num2str(i)]);
            print(gcf,'-depsc',[out_dir,out_name,'_surf_b',num2str(i)]);
        end
        % Save aggregate basis function output
        save([out_dir,out_name,'_BasisFunctions.mat'],'bvals');
    end
    % Save geometric grid data
    save([out_dir,geom.Name,'.mat'],'geom');
end
