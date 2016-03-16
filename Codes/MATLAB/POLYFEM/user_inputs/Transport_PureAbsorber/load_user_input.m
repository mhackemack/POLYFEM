function data = load_user_input(dat_in, geom_in)
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.plotSolution = 0;
data.problem.saveSolution = 0;
data.problem.saveVTKSolution = 1;
% AMR Input Parameters
% ------------------------------------------------------------------------------
data.problem.refineMesh = dat_in.refineMesh;
data.problem.refinementLevels = dat_in.refinementLevels;
data.problem.refinementTolerance = dat_in.refinementTolerance;
data.problem.AMRIrregularity = dat_in.AMRIrregularity;
data.problem.projectSolution = dat_in.projectSolution;
data.problem.refinementType = 0; % 0 = err(c)/maxerr < c, 1 = numc/totalCells = c
% Neutronics Data
% ------------------------------------------------------------------------------
data.Neutronics.PowerLevel = 1.0;
data.Neutronics.StartingSolution = 'zero';
data.Neutronics.transportMethod = 'Transport';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.SpatialMethod = dat_in.SpatialMethod;
data.Neutronics.FEMLumping = dat_in.FEMLumping;
data.Neutronics.FEMDegree = dat_in.FEMDegree;
data.Neutronics.numberEnergyGroups = 1;
% Transport Properties
% ------------------------------------------------------------------------------
% MMS Properties
data.Neutronics.Transport.MMS = true;
data.Neutronics.Transport.QuadOrder = 6;
% data.Neutronics.Transport.ExtSource = cell(data.Neutronics.numberEnergyGroups,1);
% data.Neutronics.Transport.ExactSolution = cell(data.Neutronics.numberEnergyGroups,1);
% Flux/Angle Properties
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.AngleAggregation = 'auto';
data.Neutronics.Transport.QuadType = dat_in.QuadType;
data.Neutronics.Transport.SnLevels = dat_in.SnLevels;
data.Neutronics.Transport.AzimuthalLevels = dat_in.AzimuthalLevels;
data.Neutronics.Transport.PolarLevels = dat_in.PolarLevels;
data.Neutronics.Transport.QuadAngles  = dat_in.QuadAngles;  % Angles for manual set
data.Neutronics.Transport.QuadWeights = dat_in.QuadWeights; % Weights for manual set
% Sweep Operations
data.Neutronics.Transport.performSweeps = 0;
data.Neutronics.Transport.visualizeSweeping = 0;
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Neutronics.Transport.transportType = 'upwind';
data.Neutronics.Transport.StabilizationMethod = 'EGDG';
data.Neutronics.Transport.FluxStabilization = 2.0;
data.Neutronics.Transport.CurrentStabilization = 1.0;
% Physical Properties
data.Neutronics.Transport.TotalXS = dat_in.TotalXS;
data.Neutronics.Transport.AbsorbXS = dat_in.TotalXS;
data.Neutronics.Transport.ScatteringXS = 0;
data.Neutronics.Transport.FissionXS = 0.0;
data.Neutronics.Transport.NuBar = 0.0;
data.Neutronics.Transport.FissSpec = 0.0;
data.Neutronics.Transport.ExtSource = dat_in.RHSFunc;
data.Neutronics.Transport.ExactSolution = dat_in.SolFunc;
% DSA Properties
% ------------------------------------------------------------------------------
data.Neutronics.Transport.performDSA = 0;
data.Neutronics.Transport.DSAType = 'MIP';
data.Neutronics.Transport.DSASolveMethod = 'direct';
data.Neutronics.Transport.DSAPreconditioner = 'none';
data.Neutronics.Transport.DSATolerance = 1e-4;
data.Neutronics.Transport.DSAMaxIterations = 1e7;
data.Neutronics.IP_Constant = 4;
% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.absoluteTolerance = 1e-10;
data.solver.relativeTolerance = 1e-10;
data.solver.maxIterations = 10000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];
% Get Geometry Specifications
% ------------------------------------------------------------------------------
data.problem.Dimension = geom_in.Dimension;
% [data,geometry] = load_geometry_input(data, geom_in);
% ------------------------------------------------------------------------------