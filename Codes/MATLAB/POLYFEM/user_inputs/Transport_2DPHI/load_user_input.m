function [data, geometry] = load_user_input(n,txs,c)
global glob
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.problem.Path = 'Transport/2DPHI';
data.problem.Name = 'PWLD_LS4';
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.plotSolution = 0;
data.problem.saveSolution = 0;
data.problem.saveVTKSolution = 0;
% AMR Input Parameters
% ------------------------------------------------------------------------------
data.problem.refineMesh = 0;
data.problem.refinementLevels = 8;
data.problem.refinementTolerance = 0.2;
data.problem.AMRIrregularity = 1;
data.problem.projectSolution = 0;
data.problem.refinementType = 0; % 0 = err(c)/maxerr < c, 1 = numc/totalCells = c
% Neutronics Data
% ------------------------------------------------------------------------------
data.Neutronics.PowerLevel = 1.0;
data.Neutronics.StartingSolution = 'random';
data.Neutronics.StartingSolutionFunction{1,1} = @asymptotic_limit_func;
data.Neutronics.transportMethod = 'Transport';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMLumping = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.numberEnergyGroups = 1;

% Transport Properties
% ------------------------------------------------------------------------------
% Flux/Angle Properties
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.AngleAggregation = 'auto';
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = 6;
data.Neutronics.Transport.AzimuthalLevels = 4;
data.Neutronics.Transport.PolarLevels = 2;
% Sweep Operations
data.Neutronics.Transport.performSweeps = 0;
data.Neutronics.Transport.visualizeSweeping = 0;
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Neutronics.Transport.transportType = 'upwind';
data.Neutronics.Transport.StabilizationMethod = 'EGDG';
data.Neutronics.Transport.FluxStabilization = 2.0;
data.Neutronics.Transport.CurrentStabilization = 1.0;
% Physical Properties
sigt = [txs,1/txs];
data.Neutronics.TotalXS = sigt;
data.Neutronics.DiffusionXS = (1/3)./sigt;
data.Neutronics.ScatteringXS = c*sigt;
data.Neutronics.AbsorbXS = (1-c)*sigt;
data.Neutronics.Transport.FissionXS = [0.0;0.0];
data.Neutronics.Transport.NuBar = [0.0;0.0];
data.Neutronics.Transport.FissSpec = [0.0;0.0];
data.Neutronics.Transport.ExtSource = [0.0;0.0];
% Boundary Conditions
data.Neutronics.Transport.BCFlags = [glob.Vacuum, glob.Reflecting];
data.Neutronics.Transport.BCVals  = {0.0; 0.0};

% DSA Properties
% ------------------------------------------------------------------------------
data.Neutronics.Transport.performDSA = 1;
data.Neutronics.Transport.DSAType = 'MIP';
data.Neutronics.Transport.DSASolveMethod = 'direct';
data.Neutronics.Transport.DSAPreconditioner = 'gs';
data.Neutronics.Transport.DSATolerance = 1e-3;
data.Neutronics.Transport.DSAMaxIterations = 1e3;
data.Neutronics.IP_Constant = 4;

% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.absoluteTolerance = 1e-8;
data.solver.relativeTolerance = 1e-8;
data.solver.maxIterations = 1000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 2;
L = n; ncells = n;

x=linspace(0,L,ncells+1);
y=linspace(0,L,ncells+1);
geometry = CartesianGeometry(2,x,y);

geometry.set_face_flag_on_surface(2,[0,0;0,L]);
geometry.set_face_flag_on_surface(2,[0,0;L,0]);
