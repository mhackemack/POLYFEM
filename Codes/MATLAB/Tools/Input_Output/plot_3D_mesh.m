%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Plot 3D Mesh
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
function plot_3D_mesh(mesh, cell_color, face_color, numbering)
hold on
for f=1:mesh.TotalFaces
    v = mesh.get_face_verts(f);
    vc = mesh.FaceCenter(f,:);
    if face_color
        if mesh.FaceID(f) == 0
            fill3(v(:,1), v(:,2), v(:,3),[1,1,1])
        else
            cmap = get_color_map(mesh.FaceID(f));
            fill3(v(:,1), v(:,2), v(:,3),cmap)
        end
    else
        fill3(v(:,1), v(:,2), v(:,3),[1,1,1])
    end
    if numbering == 2 || numbering == 4
        text(vc(1), vc(2), vc(3), num2str(f),'Color',[0,0,1])
    end
end
if numbering == 1 || numbering == 4
    for c=1:mesh.TotalCells
        v = mesh.CellCenter(c,:);
        text(v(1), v(2), v(3), num2str(c),'Color',[0,0,0])
    end
end
if numbering == 3 || numbering == 4
    verts = mesh.get_all_vertices();
    for V=1:mesh.TotalVertices
        v = verts(V,:);
        text(v(1), v(2), v(3), num2str(V),'Color',[1,0,0])
    end
end
if cell_color
    for c=1:mesh.TotalCells
        cmat = mesh.MatID(c);
        cmap = get_color_map(cmat);
        cfaces = mesh.CellFaces{c};
        for f=1:length(cfaces)
            ff = cfaces(f);
            fv = mesh.FaceVerts{ff};
            verts = mesh.Vertices(fv,:);
            fill3(verts(:,1), verts(:,2), verts(:,3),cmap)
        end
    end
end
hold off
alpha(.5)

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