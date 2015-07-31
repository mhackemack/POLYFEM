%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_3D_DoF(mesh, DoF)
hold on
for c=1:mesh.TotalCells
    cverts = mesh.CellVerts{c};
    cmean = mean(mesh.Vertices(cverts,:));
    text(cmean(1), cmean(2), cmean(3), num2str(c), 'Color', 'red')
end
for f=1:mesh.TotalFaces
    fcells = mesh.FaceCells(f,:);
    ffverts = mesh.Vertices(mesh.FaceVerts{f},:);
    fmean = mean(ffverts);
    fill3(ffverts(:,1),ffverts(:,2),ffverts(:,3),[1 1 1]);
    text(fmean(1), fmean(2), fmean(3), num2str(f),'Color','blue');
    if mesh.FaceID(f) ~=0, fcells(2) = []; end
    for i=1:length(fcells)
        c = fcells(i);
        fdofs = DoF.FaceCellNodes{f,i};
        fverts = DoF.NodeLocations(fdofs,:);
        cdofs = DoF.ConnectivityArray{c};
        cverts = DoF.NodeLocations(cdofs,:);
        cmean = mean(cverts);
        for j=1:length(fdofs)
            fave = (fverts(j,:) + cmean)./2;
            text(fave(1), fave(2), fave(3), num2str(fdofs(j)),'Color','black');
        end
    end
end
hold off
alpha(.5)