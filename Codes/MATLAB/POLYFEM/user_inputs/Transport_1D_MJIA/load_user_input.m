function [data, geometry] = load_user_input()
global glob
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.IO.Path = 'Transport/Homogeneous';
data.IO.Name = 'AMR_cart_Irr=3_rtol=0.2';
data.IO.PlotSolution = 0;
data.IO.SaveSolution = 0;
data.IO.SaveVTKSolution = 0;
data.IO.PrintIterationInfo = 1;
data.IO.PrintConstructionInfo = 1;
% AMR Input Parameters
% ------------------------------------------------------------------------------
data.AMR.RefineMesh = 0;
% MMS Input Parameters
% ------------------------------------------------------------------------------
data.MMS.PerformMMS = 0;
% Overall Problem Data
% ------------------------------------------------------------------------------
data.problem.NumberMaterials = 4;
data.problem.ProblemType = 'SourceDriven';
data.problem.PowerLevel = 1.0;
data.problem.TransportMethod = 'Transport';
data.problem.FEMType = 'DFEM';
data.problem.SpatialMethod = 'PWLD';
data.problem.FEMLumping = false;
data.problem.FEMDegree = 1;
% Transport Properties
% ------------------------------------------------------------------------------
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Transport.PerformAcceleration = true;
data.Transport.PnOrder = 0;
data.Transport.XSID = 1; data.Transport.QuadID = 1;
data.Transport.TransportType = 'upwind';
data.Transport.PerformSweeps = 0;
data.Transport.VisualizeSweeping = 0;
data.Transport.StabilizationMethod = 'EGDG';
data.Transport.FluxStabilization = 2.0;
data.Transport.CurrentStabilization = 1.0;
% Quadrature Properties
data.Quadrature(1).PnOrder = data.Transport.PnOrder;
data.Quadrature(1).AngleAggregation = 'auto';
data.Quadrature(1).QuadType = 'LS';
data.Quadrature(1).SnLevels = 4;
data.Quadrature(1).PolarLevels = 4;
data.Quadrature(1).AzimuthalLevels = 4;
% Flux Properties
data.Fluxes.StartingSolution = 'zero';
% Construct Group Set Information
data.Groups.NumberEnergyGroups = 99;
data.Groups.FastGroups = 1:42;
data.Groups.ThermalGroups = 43:99;
data.Groups.NumberGroupSets = length(data.Groups.FastGroups) + 1;
data.Groups.GroupSets = cell(data.Groups.NumberGroupSets,1);
data.Groups.GroupSetUpscattering = false(data.Groups.NumberGroupSets,1);
for g=1:data.Groups.NumberGroupSets, data.Groups.GroupSets{g} = g; end
data.Groups.GroupSets{end} = data.Groups.ThermalGroups;
% Retrieve All Physical Properties
% ------------------------------------------------------------------------------
% Graphite
data = add_xs_component_contribution(data, 1, 3, 'graphite_99G', 8.5238E-2);
data = add_xs_component_contribution(data, 1, 3, 'B10_99G', 2.4335449e-06);
% Air
data = add_xs_component_contribution(data, 1, 4, 'FG_CNat_99G', 7.4906E-9);
data = add_xs_component_contribution(data, 1, 4, 'N14_99G', 3.9123E-5);
data = add_xs_component_contribution(data, 1, 4, 'O16_99G', 1.0511E-5);
data = add_xs_component_contribution(data, 1, 4, 'Ar40_99G', 2.3297E-7);
% HDPE
data = add_xs_component_contribution(data, 1, 2, 'PolyH1_99G', 8.1570E-2);
data = add_xs_component_contribution(data, 1, 2, 'FG_CNat_99G', 4.0787E-2);
% BHDPE
% data = add_xs_component_contribution(data, 1, 1, 'PolyH1_99G', 5.0859E-2);
% data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 2.5429E-2);
% data = add_xs_component_contribution(data, 1, 1, 'B10_99G', 6.6256E-3);
% data = add_xs_component_contribution(data, 1, 1, 'B11_99G', 2.6669E-2);
% AmBe
data = add_xs_component_contribution(data, 1, 1, 'Am241_99G', 1.1649E-3);
data = add_xs_component_contribution(data, 1, 1, 'Be9_99G', 1.9077E-1);
data = add_xs_component_contribution(data, 1, 1, 'O16_99G', 1.0511E-5);
data.XS(1).BCFlags = [glob.Vacuum,glob.Reflecting];
data.XS(1).BCVals{1} = 0;
data.XS(1).BCVals{2} = 0;
nm = data.problem.NumberMaterials;
ng = data.Groups.NumberEnergyGroups;
data.XS(1).ExtSource = [get_PDT_AmBe_source();zeros(nm-1,ng)];
% Acceleration Properties
% ------------------------------------------------------------------------------
data.Acceleration.WGSAccelerationBool = [false(data.Groups.NumberGroupSets-1,1);1];
data.Acceleration.AGSAccelerationBool = false;
data.Acceleration.WGSAccelerationResidual = false(data.Groups.NumberGroupSets,1);
data.Acceleration.AGSAccelerationResidual = false;
data.Acceleration.WGSAccelerationID = [zeros(data.Groups.NumberGroupSets-1,1);1];
data.Acceleration.AGSAccelerationID = 0;
data.Acceleration.WGSAccelerationFrequency = ones(data.Groups.NumberGroupSets,1);
data.Acceleration.AGSAccelerationFrequency = 1;
data.Acceleration.Info(1).AccelerationType = glob.Accel_WGS_MJIA_DSA;
data.Acceleration.Info(1).DiscretizationType = glob.Accel_DSA_MIP;
data.Acceleration.Info(1).IP_Constant = 4;
data.Acceleration.Info(1).Groups = data.Groups.ThermalGroups;
data.Acceleration.Info(1).Moments = 1;
data.Acceleration.Info(1).XSID = 2;
data.Acceleration.Info(1).GroupAccelIDs = (1:length(data.Groups.ThermalGroups))+1;
data = collapse_mjia_xs(data,1,2,1);
% Loop through thermal groups and build acceleration structures
for g=1:length(data.Groups.ThermalGroups)
    tg = data.Groups.ThermalGroups(g);
    data.Acceleration.Info(g+1).AccelerationType = glob.Accel_WGS_DSA;
    data.Acceleration.Info(g+1).DiscretizationType = glob.Accel_DSA_MIP;
    data.Acceleration.Info(g+1).IP_Constant = 4;
    data.Acceleration.Info(g+1).Groups = tg;
    data.Acceleration.Info(g+1).Moments = 1;
    data.Acceleration.Info(g+1).XSID = 2+g;
    data = collapse_jacobi_xs(data,1,2+g,g+1);
end
% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.AGSMaxIterations = 1e0;
data.solver.WGSMaxIterations = 1e4*ones(data.Groups.NumberGroupSets,1);
data.solver.AGSRelativeTolerance = 1e-4;
data.solver.AGSAbsoluteTolerance = 1e-4;
data.solver.WGSRelativeTolerance = 1e-6*ones(data.Groups.NumberGroupSets,1);
data.solver.WGSAbsoluteTolerance = 1e-6*ones(data.Groups.NumberGroupSets,1);
% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 1;
L = 1e2; ncells = 40;
x=linspace(0,L,ncells+1);
geometry = CartesianGeometry(1,x);
% Set material regions
geometry.set_cell_matIDs_inside_domain(2,[5,15]);
geometry.set_cell_matIDs_inside_domain(3,[15,80]);
geometry.set_cell_matIDs_inside_domain(4,[80,100]);
% Set boundary conditions
geometry.set_face_flag_on_surface(2,0.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Listing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_PDT_AmBe_source()
out = [869835.638009
869835.638
93326.55976
88165.67049
128040.4112
39513.79114
10148.21286
163035.1443
16628.86276
107093.4906
43166.73633
267267.8753
353440.6146
37396.0072
31771.04436
13638.08021
565.4120754
131414.1515
1376.244385
1195.777371
125768.9349
49795.44883
12595.1248
2529.485137
127.48211
71.68845762
40.31326071
22.66989262
12.748211
7.168845762
4.031326071
2.266989262
1.2748211
0.716884576
0.403132607
0.226698926
0.12748211
0.071688458
0.040313261
0.022669893
0.012748211
0.007168846
0.004031326
0.002266989
0.001274821
0.000716885
0.000265914
0.000149058
0.000116298
8.85E-05
6.88E-05
3.44E-05
2.13E-05
1.44E-05
1.18E-05
9.83E-06
8.35E-06
7.37E-06
6.72E-06
5.90E-06
5.57E-06
5.08E-06
4.75E-06
4.42E-06
4.26E-06
3.93E-06
3.77E-06
3.77E-06
3.44E-06
3.12E-06
3.23E-07
3.28E-06
3.11E-06
3.11E-06
3.11E-06
2.95E-06
2.95E-06
2.54E-06
4.03E-07
2.78E-06
2.78E-06
2.95E-06
2.78E-06
2.95E-06
2.70E-07
2.68E-06
2.95E-06
3.11E-06
2.15E-07
2.90E-06
8.76E-08
3.52E-06
3.83E-06
7.28E-07
2.17E-06
1.82E-06
1.98E-06
5.16E-06
2.98E-06
]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%