function [data, geometry] = load_user_input()
global glob
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.problem.Path = 'Transport/Linear_Solution';
data.problem.Name = '';
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.plotSolution = 0;
data.problem.saveSolution = 0;
data.problem.saveVTKSolution = 1;
% AMR Input Parameters
% ------------------------------------------------------------------------------
data.problem.refineMesh = 0;
data.problem.refinementLevels = 1;
data.problem.refinementTolerance = 0.2;
data.problem.AMRIrregularity = 1;
data.problem.projectSolution = 0;
data.problem.refinementType = 0; % 0 = err(c)/maxerr < c, 1 = numc/totalCells = c
% Neutronics Data
% ------------------------------------------------------------------------------
data.Neutronics.PowerLevel = 1.0;
data.Neutronics.StartingSolution = 'zero';
data.Neutronics.transportMethod = 'Transport';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMLumping = false;
data.Neutronics.FEMDegree = 2;
data.Neutronics.numberEnergyGroups = 1;

% Transport Properties
% ------------------------------------------------------------------------------
% MMS Properties
data.Neutronics.Transport.MMS = true;
data.Neutronics.Transport.QuadOrder = 8;
data.Neutronics.Transport.ExtSource = cell(data.Neutronics.numberEnergyGroups,1);
data.Neutronics.Transport.ExactSolution = cell(data.Neutronics.numberEnergyGroups,1);
% Flux/Angle Properties
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.AngleAggregation = 'all';
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = 6;
data.Neutronics.Transport.PolarLevels = 4;
data.Neutronics.Transport.AzimuthalLevels = 4;
% Sweep Operations
data.Neutronics.Transport.performSweeps = 0;
data.Neutronics.Transport.visualizeSweeping = 0;
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Neutronics.Transport.transportType = 'upwind';
data.Neutronics.Transport.StabilizationMethod = 'EGDG';
data.Neutronics.Transport.FluxStabilization = 2.0;
data.Neutronics.Transport.CurrentStabilization = 1.0;
% Physical Properties
txs = 1.0; c = 0.0;
data.Neutronics.Transport.ScatteringXS = zeros(1,1,1,1);
data.Neutronics.Transport.TotalXS = [txs];
data.Neutronics.Transport.AbsorbXS = (1-c)*data.Neutronics.Transport.TotalXS;
data.Neutronics.Transport.ScatteringXS(1,:,:,:) = c*data.Neutronics.Transport.TotalXS;
data.Neutronics.Transport.FissionXS = [0.0];
data.Neutronics.Transport.NuBar = [0.0];
data.Neutronics.Transport.FissSpec = [0.0];
% data.Neutronics.Transport.ExtSource = [0.0];
data.Neutronics.Transport.ExtSource{1,1} = @RHSFunc_QuadSolution;
data.Neutronics.Transport.ExactSolution{1,1} = @ExactSol_QuadSolution;
% Boundary Conditions
data.Neutronics.Transport.BCFlags = [glob.Function];
data.Neutronics.Transport.BCVals{1,1} = @ExactAngSol_QuadSolution;

% DSA Properties
% ------------------------------------------------------------------------------
data.Neutronics.Transport.performDSA = 0;
data.Neutronics.Transport.DSAType = 'MIP';
data.Neutronics.IP_Constant = 4;

% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.absoluteTolerance = 1e-6;
data.solver.relativeTolerance = 1e-6;
data.solver.maxIterations = 10000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

geometry = [];
