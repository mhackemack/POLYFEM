%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Plot 1D DoF
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
function plot_1D_DoF(mesh, DoF)
hold on
xmin = mesh.minX;
xmax = mesh.maxX;
axis([xmin,xmax,-1,1])
% Loop through Vertices and plot vertical lines
for i=1:mesh.TotalVertices
    v = mesh.Vertices(i);
    plot([v,v],[-.4,.4],'k','LineWidth',2)
end
% CFEM DoF Plotter
if DoF.FEMType == 1
    for i=1:length(DoF.NodeLocations)
        v = DoF.NodeLocations(i);
        text(v,.5,num2str(i));
    end
% DFEM DoF Plotter
elseif DoF.FEMType == 2
    for c=1:mesh.TotalCells
        cn = DoF.ConnectivityArray{c};
        h = abs(DoF.NodeLocations(cn(1)) - DoF.NodeLocations(cn(2)));
        fc = mesh.CellCenter(c,:);
        % Vertex 1
        v = DoF.NodeLocations(cn(1));
        n = (fc-v) / abs(fc-v);
        text(v + n*h/4,.5,num2str(cn(1)));
        % Vertex 2
        v = DoF.NodeLocations(cn(2));
        n = (fc-v) / abs(fc-v);
        text(v + n*h/4,.5,num2str(cn(2)));
        % Remaining higher order nodes
        for i=3:length(cn)
            v = DoF.NodeLocations(cn(i));
            text(v,.5,num2str(cn(i)));
        end
    end
end
plot([xmin,xmax],[0,0],'k','LineWidth',2)
hold off
