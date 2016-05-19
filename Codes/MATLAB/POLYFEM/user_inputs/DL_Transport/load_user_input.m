function data = load_user_input()
global glob
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.problem.Path = 'Transport/IronWater';
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
data.Neutronics.StartingSolution = 'zero';
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
data.Neutronics.Transport.AngleAggregation = 'all';
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = 4;
data.Neutronics.Transport.AzimuthalLevels = 4;
data.Neutronics.Transport.PolarLevels = 2;
data.Neutronics.Transport.QuadAngles  = [1,1];  % Angles for manual set
data.Neutronics.Transport.QuadWeights = [1];  % Weights for manual set
% Sweep Operations
data.Neutronics.Transport.performSweeps = 0;
data.Neutronics.Transport.visualizeSweeping = 0;
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Neutronics.Transport.transportType = 'upwind';
data.Neutronics.Transport.StabilizationMethod = 'EGDG';
data.Neutronics.Transport.FluxStabilization = 2.0;
data.Neutronics.Transport.CurrentStabilization = 1.0;
% Physical Properties
% ep = 1e-2;
txs = 1e0; c = 0.0;
data.Neutronics.Transport.ScatteringXS = zeros(1,1,1,1);
% data.Neutronics.Transport.TotalXS = 1/ep;
% data.Neutronics.Transport.AbsorbXS = ep;
% data.Neutronics.Transport.ScatteringXS(1,:,:,:) = 1/ep-ep;
data.Neutronics.Transport.TotalXS = [txs];
data.Neutronics.Transport.AbsorbXS = (1-c)*data.Neutronics.Transport.TotalXS;
data.Neutronics.Transport.ScatteringXS(1,:,:,:) = c*data.Neutronics.Transport.TotalXS;
data.Neutronics.Transport.FissionXS = [0.0];
data.Neutronics.Transport.NuBar = [0.0];
data.Neutronics.Transport.FissSpec = [0.0];
% data.Neutronics.Transport.ExtSource = ep;
data.Neutronics.Transport.ExtSource = [1.0];
% Boundary Conditions
% data.Neutronics.Transport.BCFlags = [glob.Vacuum,glob.IncidentIsotropic];
% data.Neutronics.Transport.BCVals  = {0.0;2.0};
data.Neutronics.Transport.BCFlags = [glob.Vacuum];
data.Neutronics.Transport.BCVals  = {0.0};
% data.Neutronics.Transport.BCFlags = [glob.Vacuum, glob.IncidentIsotropic, glob.Reflecting];
% data.Neutronics.Transport.BCVals = {0.0; 1.0; 0.0};

% DSA Properties
% ------------------------------------------------------------------------------
data.Neutronics.Transport.performDSA = 0;
data.Neutronics.Transport.DSAType = 'MIP';
data.Neutronics.Transport.DSASolveMethod = 'direct';
data.Neutronics.Transport.DSAPreconditioner = 'eisenstat';
data.Neutronics.Transport.DSATolerance = 1e-4;
data.Neutronics.Transport.DSAMaxIterations = 1e3;
data.Neutronics.IP_Constant = 4;

% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.absoluteTolerance = 1e-8;
data.solver.relativeTolerance = 1e-8;
data.solver.maxIterations = 1;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = asymptotic_limit_func(x)
dim = size(x,2);
if dim == 1
    out = 0.2*cos(pi*x(:,1));
elseif dim == 2
    out = 0.2*cos(pi*x(:,1)).*cos(pi*x(:,2));
elseif dim == 3
    out = 0.2*cos(pi*x(:,1)).*cos(pi*x(:,2)).*cos(pi*x(:,3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%