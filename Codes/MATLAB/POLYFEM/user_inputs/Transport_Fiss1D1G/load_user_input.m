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
data.problem.NumberMaterials = 1;
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
data.Groups.NumberEnergyGroups = 1;
data.Groups.NumberGroupSets = data.Groups.NumberEnergyGroups;
data.Groups.GroupSets = cell(data.Groups.NumberGroupSets,1);
data.Groups.GroupSetUpscattering = false;
for g=1:data.Groups.NumberGroupSets, data.Groups.GroupSets{g} = g; end
% Retrieve All Physical Properties
data.XS(1).TotalXS      = 1.0;
data.XS(1).AbsorbXS     = 0.5;
data.XS(1).FissionXS    = 0.5;
data.XS(1).NuBar        = 2.35;
data.XS(1).FissSpec     = 1.0;
data.XS(1).ScatteringXS = 0.5;
data.XS(1).ExtSource    = 0.0;
data.XS(1).BCFlags = [glob.Vacuum];
data.XS(1).BCVals{1} = 0;
% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.AGSMaxIterations = 1;
data.solver.WGSMaxIterations = 100;
data.solver.PIMaxIterations = 1000;
data.solver.AGSRelativeTolerance = 1;
data.solver.AGSAbsoluteTolerance = 1;
data.solver.WGSRelativeTolerance = 1e-8;
data.solver.WGSAbsoluteTolerance = 1e-8;
data.solver.PIKeffTolerance = 1e-4;
data.solver.PIFluxTolerance = 1e-4;
% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 1;
L = 10; ncells = 6;
x=linspace(-L/2,L/2,ncells+1);
geometry = CartesianGeometry(1,x);

