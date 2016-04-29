function data = load_user_input()
global glob
% Problem Input Parameters
% ------------------------
data.problem.Path = 'DL_Diffusion';
data.problem.Name = '';
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.refineMesh = 0;
data.problem.refinementLevels = 0;
data.problem.refinementTolerance = 0.5;
data.problem.refinementType = 1;
data.problem.refinementSplits = 1;
data.problem.refinementType = 0; % 0 = err(c)/maxerr < c, 1 = numc/totalCells = c
data.problem.projectSolution = 0;
data.problem.plotSolution = 0;
data.problem.saveSolution = 0;

% Neutronics Data
% ---------------
data.Neutronics.PowerLevel = 1.0;
data.Neutronics.StartingSolution = 'zero';
data.Neutronics.StartingSolutionFunction{1,1} = @asymptotic_limit_func;
data.Neutronics.transportMethod = 'Diffusion';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMLumping = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.numberEnergyGroups = 1;

% Diffusion Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.Neutroncis.IP_Constant = 4;
data.Neutronics.Diffusion.MMS = false;
% Physical Properties
ep = 1e-3;
data.Neutronics.Diffusion.ScatteringXS = zeros(1,1,1);
data.Neutronics.Diffusion.DiffXS = ep/3;
data.Neutronics.Diffusion.TotalXS = 1/ep;
data.Neutronics.Diffusion.AbsorbXS = ep;
data.Neutronics.Diffusion.ScatteringXS(1,:,:) = 1/ep - ep;
% data.Neutronics.Diffusion.DiffXS = [2];
% data.Neutronics.Diffusion.TotalXS = [0];
% data.Neutronics.Diffusion.AbsorbXS = [0];
% data.Neutronics.Diffusion.ScatteringXS(1,:,:) = [0.0];
data.Neutronics.Diffusion.FissionXS = [0.0];
data.Neutronics.Transport.NuBar = [0.0];
data.Neutronics.Diffusion.FissSpec = [0.0];
data.Neutronics.Diffusion.ExtSource = ep;
% Boundary Conditions
data.Neutronics.Diffusion.BCFlags = [glob.Dirichlet];
data.Neutronics.Diffusion.BCVals = {0.0};
% data.Neutronics.Diffusion.BCFlags = [glob.Neumann; glob.Robin; glob.Robin];
% data.Neutronics.Diffusion.BCVals = {0.0;0.0;9.0};

% Solver Input Parameters
% -----------------------
data.solver.absoluteTolerance = 1e-6;
data.solver.relativeTolerance = 1e-6;
data.solver.maxIterations = 1;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Geometry Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.problem.Dimension = 2;

% Function Handles
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
