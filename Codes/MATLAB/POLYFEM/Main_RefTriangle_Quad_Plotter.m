%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Reference Polygon Quadrature Plotter
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
% Prepare Project Space
% ------------------------------------------------------------------------------
clc; close all; fclose all; format long e; clear;
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Office');
% Being User Input Section
% ------------------------------------------------------------------------------
quads = 1:8;
outdir = 'outputs/RefTriQuadPlots/';
% Execute suite and save images
% ------------------------------------------------------------------------------
verts = [0,0;1,0;0,1];
% Loop through quadrature orders
for q=1:length(quads)
    qq = quads(q);
    [qx, qw] = TriGaussPoints(qq);
    % Plot quadrature on reference triangle
    hold on; ax = gca;
    axis([0 1 0 1]); axis square;
    patch(verts(:,1), verts(:,2), [1,1,1]);
    scatter(qx(:,1), qx(:,2), 75, 'xr')
    hold off;
    ax.FontSize = 11;
    ax.XTick = linspace(0,1,6);
    ax.YTick = linspace(0,1,6);
    name = sprintf('%sRefTriQuad_Q%d',outdir,qq);
    savefig(gcf,name);
    print(gcf,'-depsc',name);
    print(gcf,'-dpng',name);
    fclose all;
end
