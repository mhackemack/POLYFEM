function [data,geometry] = load_geometry_input(data, geom_in)
data.problem.Dimension = geom_in.Dimension;
% Get Geometry
if data.problem.Dimension == 1
    geometry = CartesianGeometry(1,geom_in.x);
elseif strcmpi(geom_in.GeometryType, 'cart')
    if data.problem.Dimension == 2
        geometry = CartesianGeometry(2,geom_in.x,geom_in.y);
    elseif data.problem.Dimension == 3
        geometry = CartesianGeometry(3,geom_in.x,geom_in.y,geom_in.z);
    end
elseif strcmpi(geom_in.GeometryType, 'tri')
    [x,y]=meshgrid(geom_in.x,geom_in.y);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
    if data.problem.Dimension == 3
        
    end
elseif strcmpi(geom_in.GeometryType, 'tet')
    [x,y,z]=meshgrid(geom_in.x,geom_in.y,geom_in.z);
    x=x(:);y=y(:);z=z(:);
    tri = delaunayTriangulation(x,y,z);
    geometry = GeneralGeometry(3, 'Delaunay', tri);
elseif strcmpi(geom_in.GeometryType, 'poly')
    gname = strcat('PolyMesh_SqDomain_L1_n',num2str(geom_in.PolyNum),'.mat');
    load(gname);
elseif strcmpi(geom_in.GeometryType, 'split_poly')
    gname = strcat('SplitPolyMesh_slope=-1_yint=1_n',num2str(geom_in.PolyNum),'.mat');
    load(gname);
end
% Boundary Condition Flags
if data.problem.Dimension == 1
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type];
elseif data.problem.Dimension == 2
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type,geom_in.ymin_bound_type,geom_in.ymax_bound_type];
elseif data.problem.Dimension == 3
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type,geom_in.ymin_bound_type,geom_in.ymax_bound_type,geom_in.zmin_bound_type,geom_in.zmax_bound_type];
end
% Boundary Condition Values
data.Neutronics.Transport.BCVals = cell(2*data.problem.Dimension,1);
% xmin vals
data.Neutronics.Transport.BCVals{1} = geom_in.xmin_val;
% xmax vals
data.Neutronics.Transport.BCVals{2} = geom_in.xmax_val;
if data.problem.Dimension > 1
    % ymin vals
    data.Neutronics.Transport.BCVals{3} = geom_in.ymin_val;
    % ymax vals
    data.Neutronics.Transport.BCVals{4} = geom_in.ymax_val;
end
if data.problem.Dimension > 2
    % ymin vals
    data.Neutronics.Transport.BCVals{5} = geom_in.zmin_val;
    % ymax vals
    data.Neutronics.Transport.BCVals{6} = geom_in.zmax_val;
end
% Assign boundary flags to geometry
if data.problem.Dimension == 1
    geometry.set_face_flag_on_surface(1,0);
    geometry.set_face_flag_on_surface(2,geom_in.Lx);
elseif data.problem.Dimension == 2
    geometry.set_face_flag_on_surface(1,[0,0;0,geom_in.Ly]);
    geometry.set_face_flag_on_surface(2,[geom_in.Lx,0;geom_in.Lx,geom_in.Ly]);
    geometry.set_face_flag_on_surface(3,[0,0;geom_in.Lx,0]);
    geometry.set_face_flag_on_surface(4,[geom_in.Lx,geom_in.Ly;0,geom_in.Ly]);
elseif data.problem.Dimension == 3
    
end