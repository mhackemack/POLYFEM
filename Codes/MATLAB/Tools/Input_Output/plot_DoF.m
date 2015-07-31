function plot_DoF(mesh, DoF)
if mesh.Dimension == 1
    plot_1D_DoF(mesh, DoF)
elseif mesh.Dimension == 2
    plot_2D_DoF(mesh, DoF)
elseif mesh.Dimension == 3
    plot_3D_DoF(mesh, DoF)
end