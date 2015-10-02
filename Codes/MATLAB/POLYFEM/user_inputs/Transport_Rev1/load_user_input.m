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
% Neutronics Data
% ------------------------------------------------------------------------------
data.problem.NumberMaterials = 1;
data.problem.ProblemType = 'SourceDriven';
data.problem.PowerLevel = 1.0;
data.problem.TransportMethod = 'Transport';
data.problem.FEMType = 'DFEM';
data.problem.SpatialMethod = 'PWLD';
data.problem.FEMDegree = 1;
% Transport Properties
% ------------------------------------------------------------------------------
% Quadrature Properties
data.Quadrature(1).PnOrder = 0;
data.Quadrature(1).AngleAggregation = 'auto';
data.Quadrature(1).QuadType = 'LS';
data.Quadrature(1).SnLevels = 4;
data.Quadrature(1).PolarLevels = 4;
data.Quadrature(1).AzimuthalLevels = 4;
% Flux Properties
data.Fluxes.PnOrder = 0;
data.Fluxes.StartingSolution = 'zero';
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Transport.PerformAcceleration = true;
data.Transport.XSID = 1; data.Transport.QuadID = 1;
data.Transport.TransportType = 'upwind';
data.Transport.PerformSweeps = 0;
data.Transport.VisualizeSweeping = 0;
data.Transport.StabilizationMethod = 'EGDG';
data.Transport.FluxStabilization = 2.0;
data.Transport.CurrentStabilization = 1.0;
% Construct Group Set Information
data.Groups.NumberEnergyGroups = 1;
data.Groups.NumberGroupSets = 1;
data.Groups.GroupsSets{1} = 1;
data.Groups.GroupSetUpscattering = false;
% Retrieve All Physical Properties
txs = 1; c = 0.5;
data.XS(1).TotalXS = txs; data.XS(1).DiffXS = 1/3/txs; data.XS(1).AbsorbXS = (1-c)*txs;
data.XS(1).FissionXS = 0; data.XS(1).NuBar = 0; data.XS(1).FissSpec = 0;
data.XS(1).ScatteringXS = c*txs; data.XS(1).ExtSource = 1.0;
data.XS(1).BCFlags = [glob.Vacuum]; data.XS(1).BCVals{1} = 0;
% Acceleration Properties
% ------------------------------------------------------------------------------
data.Acceleration.WGSAccelerationBool = 0;
data.Acceleration.AGSAccelerationBool = 0;
data.Acceleration.WGSAccelerationResidual = 0;
data.Acceleration.AGSAccelerationResidual = 0;
data.Acceleration.WGSAccelerationID = 1;
data.Acceleration.AGSAccelerationID = 0;
data.Acceleration.Info(1).AccelerationType = glob.Accel_WGS_DSA;
data.Acceleration.Info(1).DiscretizationType = glob.Accel_DSA_MIP;
data.Acceleration.Info(1).IP_Constant = 4;
data.Acceleration.Info(1).Groups = 1;
data.Acceleration.Info(1).Moments = 1;
data.Acceleration.Info(1).XSID = 1;
% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.AGSMaxIterations = 1;
data.solver.WGSMaxIterations = 1000;
data.solver.AGSRelativeTolerance = 1e-8;
data.solver.WGSRelativeTolerance = 1e-8;
data.solver.AGSAbsoluteTolerance = 1e-8;
data.solver.WGSAbsoluteTolerance = 1e-8;
% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 1;
L = 1; ncells = 10;

% tx = linspace(0,L,ncells+1);
% [x,y]=meshgrid(tx,tx);
% x=x(:);y=y(:);
% tri = delaunayTriangulation(x,y);
% geometry = GeneralGeometry(2, 'Delaunay', tri);

x=linspace(0,L,ncells+1);
y=linspace(0,L,ncells+1);
z=linspace(0,L,ncells+1);
geometry = CartesianGeometry(1,x);
% geometry = CartesianGeometry(2,x,y);
% geometry = CartesianGeometry(3,x,y,z);

% geometry.extrude_mesh_2D_to_3D([0,1/3,2/3,1]);

% geometry.set_face_flag_on_surface(2,0.0);
% geometry.set_face_flag_on_surface(2,[0,.2*L;0,.4*L]);
% geometry.set_face_flag_on_surface(2,[0,L;L,L]);
% geometry.set_face_flag_on_surface(2,[L,0;L,L]);
% geometry.set_face_flag_on_surface(2,[0,0;L,0]);
