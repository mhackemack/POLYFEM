%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Plot 2D DoFHandler
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to plot the distributed 2D DoFs.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_2D_DoF(mesh, DoF)
hold on
if DoF.FEMType == 2
    for e=1:mesh.TotalCells
        vverts = mesh.Vertices(mesh.CellVerts{e},:);
        fill(vverts(:,1),vverts(:,2),[1 1 1]);
        cc = mesh.CellCenter(e,:);
        nodes = DoF.ConnectivityArray{e};
        locs = DoF.NodeLocations(nodes,:);
        for i=1:length(nodes)
            xave = (cc(1) + locs(i,1)) / 2;
            yave = (cc(2) + locs(i,2)) / 2;
            text(xave,yave,num2str(nodes(i)));
        end
    end
elseif DoF.FEMType == 1
    for e=1:mesh.TotalCells
        vverts = mesh.Vertices(mesh.CellVerts{e},:);
        fill(vverts(:,1),vverts(:,2),[1 1 1]);
    end
    for i=1:DoF.TotalDoFs
        xy = DoF.NodeLocations(i,:);
        text(xy(1), xy(2),num2str(i));
    end
end
hold off