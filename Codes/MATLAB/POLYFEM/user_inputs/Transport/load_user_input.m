function [data, geometry] = load_user_input()
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
data.Neutronics.SpatialMethod = 'LD';
data.Neutronics.FEMLumping = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.numberEnergyGroups = 1;

% Transport Properties
% ------------------------------------------------------------------------------
% Flux/Angle Properties
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.AngleAggregation = 'single';
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
data.solver.absoluteTolerance = 1e-6;
data.solver.relativeTolerance = 1e-6;
data.solver.maxIterations = 1000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 2;
L = 1; ncells = 10;
% gname = 'PolyMesh_SqDomain_L1_n256';
% gname = 'assembly_L10_4x4_R=0.6';
% gname = 'misha_quad_L1_n4';
% gname = 'random_poly_mesh_L1_n16_a0.9';
% gname = 'shestakov_poly_mesh_L1_nc4_a0.25';
% gname = 'z_mesh_quad_L1_n9_a0.15';
% gname = 'z_mesh_poly_L1_n20_a0.05';
% gname = 'smooth_quad_mesh_L1_nc5_emb6_a0.15';
% gname = 'smooth_poly_mesh_L1_n8_a0.15';
% load(strcat(glob.geom_path,gname,'.mat'));
% data = get_SimpleReactor_XS(data);

% tx = linspace(0,L,ncells+1);
% [x,y]=meshgrid(tx,tx);
% x=x(:);y=y(:);
% tri = delaunayTriangulation(x,y);
% geometry = GeneralGeometry(2, 'Delaunay', tri);

% tx = linspace(0,L,ncells+1);
% [x,y,z]=meshgrid(tx,tx,tx);
% x=x(:);y=y(:);z=z(:);
% tri = delaunayTriangulation(x,y,z);
% geometry = GeneralGeometry(3, 'Delaunay', tri);

x=linspace(0,L,ncells+1);
y=linspace(0,L,ncells+1);
% z=linspace(0,L,ncells+1);
% geometry = CartesianGeometry(1,x);
geometry = CartesianGeometry(2,x,y);
% geometry = CartesianGeometry(3,x,y,z);

% geometry.turn_2D_mesh_to_traps(.0001);

% geometry.extrude_mesh_2D_to_3D(linspace(0,L,ncells+1));
% geometry.extrude_mesh_2D_to_3D([0,1/3,2/3,1]);

% geometry.set_face_flag_on_surface(2,0.0);
% geometry.set_face_flag_on_surface(2,[0,.2*L;0,.4*L]);
% geometry.set_face_flag_on_surface(2,[0,0;0,L]);
% geometry.set_face_flag_on_surface(3,[0,L;L,L]);
% geometry.set_face_flag_on_surface(1,[L,0;L,L]);
% geometry.set_face_flag_on_surface(3,[0,0;L,0]);

% [data, geometry] = get_EIR( data, 1, 'cart' );
% [data, geometry] = get_IronWater( data, 4, 'cart' );
% [data, geometry] = get_IronWaterII( data, 1, 'tri' );
% [data, geometry] = get_BWRAssembly( data, 2, 'cart' );
% [data, geometry] = get_Yaqi_2D( data, 4, 'cart' );
% [data, geometry] = get_2D_SS_tophat( data, .9, 1, 'cart' );
% [data, geometry] = get_3D_SS_tophat( data, 1, 1, 'cart' );
% [data, geometry] = get_Reed_1D( data, 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EIR-2 Benchmark Overwrite
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
data.Neutronics.Transport.BCVals  = {0.0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iron-Water Benchmark Overwrite
function [data, geometry] = get_IronWater( data, n, s )
global glob
% Get XS
data = get_IronWater_XS( data );
% Get Geometry
if length(n) == 1, n = n*ones(1,4); end
x = 0; y = 0; nn = 5;
dx = [0, 12, 15, 21, 30];
dy = [0, 12, 15, 21, 30];
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
geometry.set_cell_matIDs_inside_domain(2, [0,0;30,0;30,30;0,30]);
geometry.set_cell_matIDs_inside_domain(3, [0,0;21,0;21,21;0,21]);
geometry.set_cell_matIDs_inside_domain(2, [0,0;15,0;15,15;0,15]);
geometry.set_cell_matIDs_inside_domain(1, [0,0;12,0;12,12;0,12]);
% Boundary Conditions
geometry.set_face_flag_on_surface(2,[0,0;0,30]);
geometry.set_face_flag_on_surface(2,[0,0;30,0]);
data.Neutronics.Transport.BCFlags = [glob.Vacuum, glob.Reflecting];
data.Neutronics.Transport.BCVals  = {0.0,         0.0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iron-Water Benchmark Overwrite
function [data, geometry] = get_IronWaterII( data, n, s )
global glob
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
% BWR Assembly Benchmark Overwrite
function [data, geometry] = get_BWRAssembly( data, n, s )
global glob
% Get XS
data = get_BWRAssembly_XS( data );
% Get Geometry
x = 0; y = 0; nn = 12;
dx = [.47498,.34544,1.87452,1.87452,1.87452,1.87452,1.87452,1.87452,1.87452,.34544,.47625,.47625];
dy = [.47498,.34544,1.87452,1.87452,1.87452,1.87452,1.87452,1.87452,1.87452,.34544,.47625,.47625];
xx = 0; yy = 0; xm = 0; ym = 0;
for i=1:nn
    xt = linspace(xx, xx+dx(i), n+1);
    yt = linspace(yy, yy+dy(i), n+1);
    x = [x, xt(2:end)];
    y = [y, yt(2:end)];
    xx = xx + dx(i); yy = yy + dy(i);
    xm(i+1) = xm(i) + dx(i);
    ym(i+1) = ym(i) + dy(i);
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
geometry.set_cell_matIDs_inside_domain(7, [0,0;xm(end),0;xm(end),ym(end);0,ym(end)]);
geometry.set_cell_matIDs_inside_domain(6, [xm(2),ym(2);xm(end-2),ym(2);xm(end-2),ym(end-2);xm(2),ym(end-2)]);
geometry.set_cell_matIDs_inside_domain(1, [xm(3),ym(3);xm(end-3),ym(3);xm(end-3),ym(end-3);xm(3),ym(end-3)]);
geometry.set_cell_matIDs_inside_domain(3, [xm(3),ym(end-4);xm(end-4),ym(end-4);xm(end-4),ym(end-3);xm(3),ym(end-3)]);
geometry.set_cell_matIDs_inside_domain(3, [xm(end-4),ym(3);xm(end-3),ym(3);xm(end-3),ym(end-4);xm(end-4),ym(end-4)]);
geometry.set_cell_matIDs_inside_domain(4, [xm(end-4),ym(end-4);xm(end-3),ym(end-4);xm(end-3),ym(end-3);xm(end-4),ym(end-3)]);
geometry.set_cell_matIDs_inside_domain(2, [xm(3),ym(3);xm(4),ym(3);xm(4),ym(4);xm(3),ym(4)]);
geometry.set_cell_matIDs_inside_domain(2, [xm(4),ym(end-4);xm(7),ym(end-4);xm(7),ym(end-3);xm(4),ym(end-3)]);
geometry.set_cell_matIDs_inside_domain(2, [xm(end-4),ym(4);xm(end-3),ym(4);xm(end-3),ym(7);xm(end-4),ym(7)]);
geometry.set_cell_matIDs_inside_domain(2, [xm(end-5),ym(end-5);xm(end-4),ym(end-5);xm(end-4),ym(end-4);xm(end-5),ym(end-4)]);
geometry.set_cell_matIDs_inside_domain(5, [xm(6),ym(4);xm(7),ym(4);xm(7),ym(5);xm(6),ym(5)]);
geometry.set_cell_matIDs_inside_domain(5, [xm(4),ym(6);xm(5),ym(6);xm(5),ym(7);xm(4),ym(7)]);
geometry.set_cell_matIDs_inside_domain(5, [xm(8),ym(6);xm(9),ym(6);xm(9),ym(7);xm(8),ym(7)]);
geometry.set_cell_matIDs_inside_domain(5, [xm(6),ym(8);xm(7),ym(8);xm(7),ym(9);xm(6),ym(9)]);
% Boundary Conditions
geometry.set_face_flag_on_surface(2,[0,0;0,y(end)]);
geometry.set_face_flag_on_surface(2,[0,0;x(end),0]);
data.Neutronics.Transport.BCFlags = [glob.Vacuum, glob.Reflecting];
data.Neutronics.Transport.BCVals  = {0.0,         0.0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, geometry] = get_Yaqi_2D( data, n, s )
global glob
% Get XS
data = get_Yaqi_XS( data );
% Get Geometry
ntot = n*5 + 1;
Lx = 2.0; Ly = 1.0;
dx = Lx / (5); dy = Ly / (5);
if strcmp(s, 'cart')
    geometry = CartesianGeometry(2,linspace(0,Lx,ntot),linspace(0,Ly,ntot));
elseif strcmp(s, 'tri')
    [x,y]=meshgrid(linspace(0,Lx,ntot),linspace(0,Ly,ntot));
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
end
% Set Material IDs
geometry.set_cell_matIDs_inside_domain(4, [2*dx,2*dy;Lx,2*dy;Lx,Ly;2*dx,Ly]);
geometry.set_cell_matIDs_inside_domain(1, [0,0;dx,0;dx,Ly;0,Ly]);
geometry.set_cell_matIDs_inside_domain(2, [dx,0;3*dx,0;3*dx,2*dy;dx,2*dy]);
geometry.set_cell_matIDs_inside_domain(3, [3*dx,0;Lx,0;Lx,2*dy;3*dx,2*dy]);
geometry.set_cell_matIDs_inside_domain(3, [dx,2*dy;2*dx,2*dy;2*dx,Ly;dx,Ly]);
% Boundary Conditions
data.Neutronics.Transport.BCFlags = glob.Vacuum;
data.Neutronics.Transport.BCVals  = {0.0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, geometry] = get_2D_SS_tophat( data, c, n, s )
global glob
data.problem.Dimension = 2;
% Get XS
data = get_SS_TopHat_XS( data, c );
% Get Geometry
Lx = 7; Ly = 4;
nx = n*14; ny = n*8;
x = linspace(0,Lx,nx+1); y = linspace(0,Ly,ny+1);
if strcmp(s, 'cart')
    geometry = CartesianGeometry(2,x,y);
elseif strcmp(s, 'tri')
    [x,y]=meshgrid(x,y);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
end
% Set Material IDs
geometry.set_cell_matIDs_inside_domain(1, [0,0;Lx,0;Lx,Ly;0,Ly]);
geometry.set_cell_matIDs_inside_domain(2, [0,1.5;Lx,1.5;Lx,2.5;0,2.5]);
geometry.set_cell_matIDs_inside_domain(2, [2.5,0.5;4.5,0.5;4.5,3.5;2.5,3.5]);
geometry.set_cell_matIDs_inside_domain(1, [3.0,1.0;4.0,1.0;4.0,3.0;3.0,3.0]);
% Boundary Conditions
geometry.set_face_flag_on_surface(2,[0,1.5;0,2.5]);
data.Neutronics.Transport.BCFlags = [glob.Vacuum, glob.IncidentIsotropic];
data.Neutronics.Transport.BCVals  = {0.0, 1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, geometry] = get_3D_SS_tophat( data, c, n, s )
global glob
data.problem.Dimension = 3;
% Get XS
data = get_SS_TopHat_XS( data, c );
% Get Geometry
Lx = 7; Ly = 4; Lz = 4;
nx = n*14; ny = n*8; nz = n*8;
x = linspace(0,Lx,nx+1); y = linspace(0,Ly,ny+1); z = linspace(0,Lz,nz+1);
if strcmp(s, 'cart')
    geometry = CartesianGeometry(2,x,y);
elseif strcmp(s, 'tri')
    [x,y]=meshgrid(x,y);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
end
% Set Material IDs
geometry.set_cell_matIDs_inside_domain(1, [0,0;Lx,0;Lx,Ly;0,Ly]);
geometry.set_cell_matIDs_inside_domain(2, [0,1.5;Lx,1.5;Lx,2.5;0,2.5]);
geometry.set_cell_matIDs_inside_domain(2, [2.5,0.5;4.5,0.5;4.5,3.5;2.5,3.5]);
geometry.set_cell_matIDs_inside_domain(1, [3.0,1.0;4.0,1.0;4.0,3.0;3.0,3.0]);
% Boundary Conditions
geometry.set_face_flag_on_surface(2,[0,1.5;0,2.5]);
data.Neutronics.Transport.BCFlags = [glob.Vacuum, glob.IncidentIsotropic];
data.Neutronics.Transport.BCVals  = {0.0, 1e5};
% Extrude Mesh
geometry.extrude_mesh_2D_to_3D(z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, geometry] = get_Reed_1D( data, n )
global glob
data.problem.Dimension = 1;
% Get XS
data = get_Reed_XS( data );
% Get Geometry
L = 8; nx = n*8;
x = linspace(0,L,nx+1);
geometry = CartesianGeometry(1,x);
% Set Material IDs
geometry.set_cell_matIDs_inside_domain(1, [0,2]);
geometry.set_cell_matIDs_inside_domain(2, [2,3]);
geometry.set_cell_matIDs_inside_domain(3, [3,5]);
geometry.set_cell_matIDs_inside_domain(4, [5,6]);
geometry.set_cell_matIDs_inside_domain(5, [6,8]);
% Boundary Conditions
data.Neutronics.Transport.BCFlags = [glob.Vacuum];
data.Neutronics.Transport.BCVals  = {0.0};
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