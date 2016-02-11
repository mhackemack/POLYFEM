function [data, geometry] = load_user_input(dat_in)
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.plotSolution = 0;
data.problem.saveSolution = 1;
data.problem.saveVTKSolution = 1;
% AMR Input Parameters
% ------------------------------------------------------------------------------
data.problem.refineMesh = 1;
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
% Flux/Angle Properties
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.AngleAggregation = 'single';
data.Neutronics.Transport.QuadType = dat_in.QuadType;
data.Neutronics.Transport.SnLevels = dat_in.SnLevels;
data.Neutronics.Transport.AzimuthalLevels = dat_in.AzimuthalLevels;
data.Neutronics.Transport.PolarLevels = dat_in.PolarLevels;
% Sweep Operations
data.Neutronics.Transport.performSweeps = 0;
data.Neutronics.Transport.visualizeSweeping = 0;
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Neutronics.Transport.transportType = 'upwind';
data.Neutronics.Transport.StabilizationMethod = 'EGDG';
data.Neutronics.Transport.FluxStabilization = 2.0;
data.Neutronics.Transport.CurrentStabilization = 1.0;
% DSA Properties
% ------------------------------------------------------------------------------
data.Neutronics.Transport.performDSA = 1;
data.Neutronics.Transport.DSAType = 'MIP';
data.Neutronics.Transport.DSASolveMethod = dat_in.DSASolveMethod;
data.Neutronics.Transport.DSAPreconditioner = dat_in.DSAPreconditioner;
data.Neutronics.Transport.DSATolerance = dat_in.DSATolerance;
data.Neutronics.Transport.DSAMaxIterations = 1e7;
data.Neutronics.IP_Constant = 4;
% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.absoluteTolerance = 1e-8;
data.solver.relativeTolerance = 1e-8;
data.solver.maxIterations = 10000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];
% Get IronWater Specifications
% ------------------------------------------------------------------------------
[data, geometry] = get_IronWaterII( data, 1, dat_in.GeometryType );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iron-Water Benchmark Overwrite
function [data, geometry] = get_IronWaterII( data, n, s )
global glob
data.problem.Dimension = 2;
% Get XS
data = get_IronWaterII_XS( data );
% Get Geometry
L = 10;
x = linspace(0,L,5*n+1);
if strcmp(s, 'cart')
    geometry = CartesianGeometry(2,x,x);
elseif strcmp(s, 'tri')
    [x,y]=meshgrid(x,x);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
end
% Set Material IDs
geometry.set_cell_matIDs_inside_domain(3, [0,0;L,0;L,L;0,L]);
geometry.set_cell_matIDs_inside_domain(2, [0,6;4,6;4,10;0,10]);
geometry.set_cell_matIDs_inside_domain(1, [0,8;2,8;2,10;0,10]);
% Boundary Conditions
geometry.set_face_flag_on_surface(2,[0,0;0,L]);
geometry.set_face_flag_on_surface(2,[0,L;L,L]);
data.Neutronics.Transport.BCFlags = [glob.Vacuum, glob.Reflecting];
data.Neutronics.Transport.BCVals  = {0.0,         0.0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%