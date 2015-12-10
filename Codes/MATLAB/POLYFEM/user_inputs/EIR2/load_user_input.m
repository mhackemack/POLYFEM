function [data, geometry] = load_user_input(dat_in)
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.plotSolution = 0;
data.problem.saveSolution = 0;
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
data.Neutronics.FEMDegree = dat_in.FEMDegree;
data.Neutronics.numberEnergyGroups = 1;
% Transport Properties
% ------------------------------------------------------------------------------
% Flux/Angle Properties
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.AngleAggregation = 'auto';
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
data.Neutronics.Transport.DSAPreconditioner = 'gs';
data.Neutronics.Transport.DSATolerance = 1e-4;
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
[data, geometry] = get_EIR( data, 1, dat_in.GeometryType );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iron-Water Benchmark Overwrite
function [data, geometry] = get_EIR( data, n, s )
global glob
% Get XS
data = get_EIR2_XS( data );
% Get Geometry
data.problem.Dimension = 2;
if length(n) == 1, n = n*ones(1,4); end
x = 0; y = 0; nn = 5;
dx = [0, 18, 48, 78, 96];
dy = [0, 18, 43, 68, 86];
for i=1:nn-1
    xt = linspace(dx(i), dx(i+1), n(i)+1);
    yt = linspace(dy(i), dy(i+1), n(i)+1);
    x = [x, xt(2:end)];
    y = [y, yt(2:end)];
end
if strcmp(s, 'cart')
    geometry = CartesianGeometry(2,x,y);
elseif strcmp(s, 'tri')
    [x,y]=meshgrid(x,y);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
end
% Set Material IDs
geometry.set_cell_matIDs_inside_domain(2, [18,18;48,18;48,43;18,43]);
geometry.set_cell_matIDs_inside_domain(3, [48,18;78,18;78,43;48,43]);
geometry.set_cell_matIDs_inside_domain(4, [48,43;78,43;78,68;48,68]);
geometry.set_cell_matIDs_inside_domain(5, [18,43;48,43;48,68;18,68]);
% Boundary Conditions
data.Neutronics.Transport.BCFlags = glob.Vacuum;
data.Neutronics.Transport.BCVals  = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%