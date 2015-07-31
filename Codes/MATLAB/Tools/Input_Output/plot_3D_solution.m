function plot_3D_solution(mesh, DoF, ~, x, basis_name)

if iscell(x), x = x{1}; end
mtype = mesh.get_mesh_type();
ftype = DoF.FEMName;

if strcmp(ftype, 'CFEM')
    if strcmp(mtype, 'Tetrahedron')
        plot_tetrahedron_cfem_solution(mesh,DoF,x);
    else
        plot_polyhedral_cfem_solution(mesh,DoF,x);
    end
else
    plot_polyhedral_dfem_solution(mesh,DoF,x);
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tetrahedron_cfem_solution(mesh,DoF,sol)
verts = mesh.Vertices;
x = verts(:,1);
y = verts(:,2);
z = verts(:,3);
surf(x,y,z,sol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_polyhedral_cfem_solution(mesh,DoF,sol)
for f=1:mesh.TotalFaces
    fflag = mesh.get_face_flags(f);
    v = mesh.get_face_verts(f);
    if fflag == 0
        d1 = DoF.FaceCellNodes{f,1};
        d2 = DoF.FaceCellNodes{f,2};
        patch(v(:,1), v(:,2), v(:,3), sol(d1))
        patch(v(:,1), v(:,2), v(:,3), sol(d2))
    else
        d = DoF.FaceCellNodes{f,1};
        patch(v(:,1), v(:,2), v(:,3), sol(d))
    end
end
alpha(.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tetrahedron_dfem_solution(mesh,DoF,sol)

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_polyhedral_dfem_solution(mesh,DoF,sol)
for f=1:mesh.TotalFaces
    fflag = mesh.get_face_flags(f);
    v = mesh.get_face_verts(f);
    if fflag == 0
        d1 = DoF.FaceCellNodes{f,1};
        d2 = DoF.FaceCellNodes{f,2};
        patch(v(:,1), v(:,2), v(:,3), sol(d1))
        patch(v(:,1), v(:,2), v(:,3), sol(d2))
    else
        d = DoF.FaceCellNodes{f,1};
        patch(v(:,1), v(:,2), v(:,3), sol(d))
    end
end
alpha(.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
