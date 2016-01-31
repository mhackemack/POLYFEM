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
data.problem.NumberMaterials = 2;
data.problem.ProblemType = 'Eigenvalue';
data.problem.KeffGuess = 1.0;
data.problem.PowerLevel = 1.0;
data.problem.TransportMethod = 'Transport';
data.problem.FEMType = 'DFEM';
data.problem.SpatialMethod = 'PWLD';
data.problem.FEMLumping = false;
data.problem.FEMDegree = 1;
% Transport Properties
% ------------------------------------------------------------------------------
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Transport.PerformAcceleration = false;
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
data.Quadrature(1).AngleAggregation = 'all';
data.Quadrature(1).QuadType = 'LS';
data.Quadrature(1).SnLevels = 8;
% Flux Properties
data.Fluxes.StartingSolution = 'one';
% Construct Group Set Information
data.Groups.NumberEnergyGroups = 7;
data.Groups.FastGroups = 1:3;
data.Groups.ThermalGroups = 4:7;
data.Groups.NumberGroupSets = 1;
data.Groups.GroupSets = cell(data.Groups.NumberGroupSets,1);
data.Groups.GroupSetUpscattering = false;
data.Groups.GroupSets{1} = 1:data.Groups.NumberEnergyGroups;
% data.Groups.NumberGroupSets = data.Groups.NumberEnergyGroups;
% data.Groups.GroupSets = cell(data.Groups.NumberGroupSets,1);
% data.Groups.GroupSetUpscattering = [false(length(data.Groups.FastGroups),1);true(length(data.Groups.ThermalGroups),1)];
% for g=1:data.Groups.NumberGroupSets, data.Groups.GroupSets{g} = g; end
% Retrieve All Physical Properties
% data.XS(1) = get_C5G7_XS(1);
data.XS(1) = get_C5G7_XS([7,1]);
data.XS(1).BCFlags = [glob.Vacuum,glob.Reflecting];
data.XS(1).BCVals{1} = 0;
% Acceleration Properties
% ------------------------------------------------------------------------------
data.Acceleration.WGSAccelerationBool = false;
data.Acceleration.AGSAccelerationBool = false;
data.Acceleration.WGSAccelerationResidual = false(data.Groups.NumberGroupSets,1);
data.Acceleration.AGSAccelerationResidual = false;
data.Acceleration.WGSAccelerationID = 1;
data.Acceleration.AGSAccelerationID = 0;
data.Acceleration.WGSAccelerationFrequency = ones(data.Groups.NumberGroupSets,1);
data.Acceleration.AGSAccelerationFrequency = 1;
data.Acceleration.Info(1).AccelerationType = glob.Accel_WGS_DSA;
data.Acceleration.Info(1).DiscretizationType = glob.Accel_DSA_MIP;
data.Acceleration.Info(1).IP_Constant = 4;
data.Acceleration.Info(1).Groups = 1:data.Groups.NumberEnergyGroups;
data.Acceleration.Info(1).Moments = 1;
data.Acceleration.Info(1).XSID = 2;
data = collapse_jacobi_xs(data,1,2,1);
% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.AGSMaxIterations = 1;
data.solver.WGSMaxIterations = 10*ones(data.Groups.NumberGroupSets,1);
data.solver.PIMaxIterations = 10000;
data.solver.AGSRelativeTolerance = 1e-1;
data.solver.AGSAbsoluteTolerance = 1e-1;
data.solver.WGSRelativeTolerance = 1e-3*ones(data.Groups.NumberGroupSets,1);
data.solver.WGSAbsoluteTolerance = 1e-3*ones(data.Groups.NumberGroupSets,1);
data.solver.PIKeffTolerance = 1e-5;
data.solver.PIFluxTolerance = 1e-5;
% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 2;
L = 10; ncells = 4;
x=linspace(0,L,ncells+1);
y=linspace(0,L,ncells+1);
geometry = CartesianGeometry(2,x,y);
% Set material regions
geometry.set_cell_matIDs_inside_domain(2,[2.5,2.5;7.5,2.5;7.5,7.5;2.5,7.5]);
% Set boundary conditions
geometry.set_face_flag_on_surface(2,[0,0;0,10]);
geometry.set_face_flag_on_surface(2,[10,10;0,10]);