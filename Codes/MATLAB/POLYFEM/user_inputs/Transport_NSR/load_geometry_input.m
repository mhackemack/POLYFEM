function geometry = load_geometry_input(dim, m_type, dx, num)
global glob
if dim == 1
    geometry = CartesianGeometry(1,dx);
elseif dim == 2 && strcmp(m_type, 'tri')
    [x,y]=meshgrid(dx,dx);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
elseif dim == 2 && strcmp(m_type, 'quad')
    geometry = CartesianGeometry(2,dx,dx);
elseif dim == 3 && strcmp(m_type, 'hex')
    geometry = CartesianGeometry(3,dx,dx,dx);
elseif dim == 3 && strcmp(m_type, 'tet')
    [x,y,z]=meshgrid(dx,dx,dx);
    x=x(:);y=y(:);z=z(:);
    tri = delaunayTriangulation(x,y,z);
    geometry = GeneralGeometry(3, 'Delaunay', tri);
elseif dim == 3 && strcmp(m_type, 'tri_prism')
    [x,y]=meshgrid(dx,dx);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
    geometry.extrude_mesh_2D_to_3D(dx);
elseif dim == 2 && strcmp(m_type, 'poly')
    gname = ['random_poly_mesh_L1_n', num2str(2^num), '_a0.9'];
    load(strcat(glob.geom_path, gname, '.mat'));
    geometry.calculate_orthogonal_projections();
elseif dim == 3 && strcmp(m_type, 'poly')
    L = max(dx);
    gname = ['random_poly_mesh_L1_n', num2str(2^num), '_a0.9'];
    load(strcat(glob.geom_path, gname, '.mat'));
    dxy = mean(geometry.CellVolume)^(1/2);
    nz = floor(L/dxy)+1;
    z = linspace(0,L,nz);
    geometry.extrude_mesh_2D_to_3D(z);
else
    error('Could not determine geometry type.')
end