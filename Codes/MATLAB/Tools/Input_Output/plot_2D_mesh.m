%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Plot 2D Mesh
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
function plot_2D_mesh(mesh,coloring,naming)
hold on
if nargin < 2
    coloring = 0;
end
if nargin < 3
    naming = 0;
end
% Loop through Cells
for e=1:mesh.TotalCells
    vverts = mesh.Vertices(mesh.CellVerts{e},:);
    if coloring == 1
        fill(vverts(:,1),vverts(:,2),get_color_map(mesh.MatID(e)));
    else
        fill(vverts(:,1),vverts(:,2),[1 1 1]);
    end
    if naming == 1 || naming == 4
        cf = mesh.CellCenter(e,:);
        text(cf(1),cf(2),num2str(e),'Color',[0,0,0]);
    end
end
% Loop through Edges
for e=1:mesh.TotalFaces
    vverts = mesh.Vertices(mesh.FaceVerts{e},:);
    plot(vverts(:,1),vverts(:,2),'k');
    if naming == 2 || naming == 4
        cf = mesh.FaceCenter(e,:);
        fflag = mesh.FaceID(e);
        text(cf(1),cf(2),num2str(e),'Color',[0,0,1]);
    end
end
% Loop through Vertices
for e=1:mesh.TotalVertices
    if naming == 3 || naming == 4
        cf = mesh.Vertices(e,:);
        text(cf(1),cf(2),num2str(e),'Color',[1,0,0]);
    end
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_color_map(val)

switch(val)
    case (1)
        out = [0 0 1];
    case (2)
        out = [0 1 0];
    case (3)
        out = [1 0 0];
    case (4)
        out = [1 0 1];
    case (5)
        out = [0 1 1];
    case (6)
        out = [1 1 0];
    case (7)
        out = [0.5 0.2 0.6];
    otherwise
        out = [0 0 1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%