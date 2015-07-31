%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Plot 1D Mesh
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
function plot_1D_mesh(mesh, coloring, naming)
hold on
xmin = mesh.minX;
xmax = mesh.maxX;
axis([xmin,xmax,-1,1])
if nargin < 2
    coloring = 0;
end
if nargin < 3
    naming = 0;
end
% Loop through Vertices and plot vertical lines
for i=1:mesh.TotalVertices
    v = mesh.Vertices(i);
    plot([v,v],[-.4,.4],'k','LineWidth',2)
    if naming == 3 || naming == 4
        text(v,.5,num2str(i),'Color',[1,0,0]);
    end
end
% Loop through faces and write text if necessary
if naming == 2 || naming == 4
    for f=1:mesh.TotalFaces
        fv = mesh.FaceCenter(f);
        text(fv,-.5,num2str(f),'Color',[0,0,1]);
    end
end
% Loop through cells
for c=1:mesh.TotalCells
    vv = mesh.CellVerts{c};
    verts = mesh.Vertices(vv,:);
    mID = mesh.MatID(c);
    if coloring
        
    else
        plot(verts,[0,0],'ko-','MarkerSize',8)
    end
    if naming == 1 || naming == 4
        cc = mesh.CellCenter(c);
        text(cc,.2,num2str(c),'Color',[0,0,0]);
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