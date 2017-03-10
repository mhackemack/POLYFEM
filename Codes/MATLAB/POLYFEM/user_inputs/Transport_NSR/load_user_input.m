function data = load_user_input(dim, bf, fdeg, bc_type)
global glob
% Problem Input Parameters
% ------------------------
data.problem.Dimension = dim;
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.refineMesh = 0;
data.problem.refinementLevels = 0;
data.problem.refinementTolerance = 0.5;
data.problem.refinementType = 1;
data.problem.refinementSplits = 1;
data.problem.plotSolution = 0;

% Neutronics Data
% ---------------
data.Neutronics.transportMethod = 'Transport';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.SpatialMethod = bf;
data.Neutronics.FEMDegree = fdeg;
data.Neutronics.FEMLumping = false;
data.Neutronics.numberEnergyGroups = 1;
data.Neutronics.StartingSolution = 'random';

% Transport Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transport Properties

% Boundary Conditions
data.Neutronics.Transport.BCVals  = {0.0};
if strcmp(bc_type, 'Vacuum')
    data.Neutronics.Transport.BCFlags = [glob.Vacuum];
elseif strcmp(bc_type, 'Reflecting')
    data.Neutronics.Transport.BCFlags = [glob.Reflecting];
end
% DSA Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.Neutronics.Transport.performDSA = 1;
data.Neutronics.Transport.DSAType = 'MIP';
data.Neutronics.Transport.DSASolveMethod = 'direct';
data.Neutronics.Transport.DSAPreconditioner = 'Jacobi';
data.Neutronics.Transport.DSATolerance = 1e-4;
data.Neutronics.Transport.DSAMaxIterations = 1e4;
data.Neutronics.IP_Constant = 4;

% Solver Input Parameters
% -----------------------
data.solver.absoluteTolerance = 1e-80;
data.solver.relativeTolerance = 1e-80;
data.solver.maxIterations = 20;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];
