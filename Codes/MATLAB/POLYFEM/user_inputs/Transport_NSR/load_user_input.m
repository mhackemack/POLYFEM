function data = load_user_input(dim, bf, bc_type)
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
data.Neutronics.FEMDegree = 1;
data.Neutronics.numberEnergyGroups = 1;

% Transport Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transport Properties

% Boundary Conditions
data.Neutronics.Transport.BCVals  = [0.0];
if strcmp(bc_type, 'Vacuum')
    data.Neutronics.Transport.BCFlags = [glob.Vacuum];
elseif strcmp(bc_type, 'Reflecting')
    data.Neutronics.Transport.BCFlags = [glob.Reflecting];
end
% DSA Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.Neutronics.Transport.performDSA = 1;
data.Neutronics.Transport.DSAType = 'MIP';
data.Neutronics.IP_Constant = 4;

% Solver Input Parameters
% -----------------------
data.solver.absoluteTolerance = 1e-10;
data.solver.relativeTolerance = 1e-8;
data.solver.maxIterations = 24;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];
