function out = get_path()

out = [];
out = [out, genpath('geometry_inputs')];
out = [out, genpath('../Tools')];
out = [out, genpath('../Extensive_Tools/export_fig')];
out = [out, genpath('../Extensive_Tools/geom3d')];
out = [out, genpath('../Extensive_Tools/geom2d')];
out = [out, genpath('../Extensive_Tools/max_ent')];
out = [out, genpath('../Extensive_Tools/Mesh2d')];
% out = [out, genpath('../Extensive_Tools/mexcdf')];
out = [out, genpath('../Extensive_Tools/PolyMesher')];
out = [out, genpath('../Extensive_Tools/Quadrature')];
out = [out, genpath('../Extensive_Tools/viz_random')];
% Add Java Files
% javaaddpath('../Extensive_Tools/javafiles/netcdfAll-4.3.jar');