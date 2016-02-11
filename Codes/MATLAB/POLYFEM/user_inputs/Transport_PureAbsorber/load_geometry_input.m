function [data,geometry] = load_geometry_input(data, geom_in)
global glob
% Get Geometry
if data.problem.Dimension == 1
    x = linspace(0,geom_in.Lx,geom_in.ncellx+1);
    geometry = CartesianGeometry(1,x);
elseif strcmpi(geom_in.GeometryType, 'cart')
    x = linspace(0,geom_in.Lx,geom_in.ncellx+1);
    y = linspace(0,geom_in.Ly,geom_in.ncelly+1);
    if data.problem.Dimension == 2
        geometry = CartesianGeometry(2,x,y);
    elseif data.problem.Dimension == 3
        z = linspace(0,geom_in.Lz,geom_in.ncellz+1);
        geometry = CartesianGeometry(3,x,y,z);
    end
elseif strcmpi(geom_in.GeometryType, 'tri')
    x = linspace(0,geom_in.Lx,geom_in.ncellx+1);
    y = linspace(0,geom_in.Ly,geom_in.ncelly+1);
    [x,y]=meshgrid(x,y);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
    if data.problem.Dimension == 3
        
    end
elseif strcmpi(geom_in.GeometryType, 'tet')
    x = linspace(0,geom_in.Lx,geom_in.ncellx+1);
    y = linspace(0,geom_in.Ly,geom_in.ncelly+1);
    z = linspace(0,geom_in.Lz,geom_in.ncellz+1);
    [x,y,z]=meshgrid(x,y,z);
    x=x(:);y=y(:);z=z(:);
    tri = delaunayTriangulation(x,y,z);
    geometry = GeneralGeometry(3, 'Delaunay', tri);
elseif strcmpi(geom_in.GeometryType, 'poly')
    gname = ['PolyMesh_SqDomain_L1_n',num2str(geom_in.PolyNum),'.mat'];
end
% Boundary Conditions
if data.problem.Dimension == 1
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type];
    data.Neutronics.Transport.BCVals = [0.0];
elseif data.problem.Dimension == 2
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type,geom_in.ymin_bound_type,geom_in.ymax_bound_type];
    
elseif data.problem.Dimension == 3
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type,geom_in.ymin_bound_type,geom_in.ymax_bound_type,geom_in.zmin_bound_type,geom_in.zmax_bound_type];
    
end