function out = get_path()

out = [];
out = [out, genpath('../Tools')];
out = [out, genpath('../Extensive_Tools/export_fig')];
out = [out, genpath('../Extensive_Tools/geom_3d')];
out = [out, genpath('../Extensive_Tools/geom_2d')];
out = [out, genpath('../Extensive_Tools/Mesh2d')];
out = [out, genpath('../Extensive_Tools/PolyMesher')];
out = [out, genpath('../Extensive_Tools/viz_random')];