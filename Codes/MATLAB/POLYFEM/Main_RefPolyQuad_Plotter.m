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
polynums = 4;
quads = 1:8;
outdir = 'outputs/RefPolyQuadPlots/';
% Execute suite and save images
% ------------------------------------------------------------------------------
% Loop through polygons
for p=1:length(polynums)
    nv = polynums(p);
    [gv,gf] = RegularPolygon(nv,1); gc = mean(gv);
    % Loop through quadrature orders
    for q=1:length(quads)
        qq = quads(q);
        [qx, qw] = get_general_volume_quadrature(gv, gf, qq, true);
        % Plot quadrature on reference polygon
        hold on;
        axis([-1 1 -1 1]); axis square;
        patch(gv(:,1),gv(:,2),[1,1,1]);
        scatter(qx(:,1), qx(:,2), 'xr')
        % Loop through vertices and apply dashed lines
        for i=1:nv
            plot([gv(i,1);gc(1)],[gv(i,2);gc(2)],'--k')
        end
        hold off;
        % Save output plots
        name = sprintf('%sV%d_Q%d',outdir,nv,qq);
        savefig(gcf,name);
        print(gcf,'-depsc',name);
        print(gcf,'-dpng',name);
    end
    fclose all;
end