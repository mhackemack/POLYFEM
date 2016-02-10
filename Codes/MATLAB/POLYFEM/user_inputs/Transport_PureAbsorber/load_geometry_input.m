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
    [x,y]=meshgrid(linspace(0,1,n+1));
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
elseif strcmpi(geom_in.GeometryType, 'poly')
    
end
% Boundary Conditions
data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type];
data.Neutronics.Transport.BCVals = [0.0,];
if data.problem.Dimension == 1
    
elseif data.problem.Dimension == 2
    
elseif data.problem.Dimension == 3
    
end