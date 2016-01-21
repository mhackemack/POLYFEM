function [data, geometry] = load_user_input()
global glob
% Problem Input Parameters
% ------------------------
data.problem.Path = 'Diffusion';
data.problem.Name = 'Diffusion_Testing';
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
data.Neutronics.StartingSolution = 'function';
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
data.Neutronics.Diffusion.BCVals = [0.0];
% data.Neutronics.Diffusion.BCFlags = [glob.Neumann; glob.Robin; glob.Robin];
% data.Neutronics.Diffusion.BCVals = [0.0;0.0;9.0];

% Solver Input Parameters
% -----------------------
data.solver.absoluteTolerance = 1e-6;
data.solver.relativeTolerance = 1e-6;
data.solver.maxIterations = 50000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Geometry Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1; n=20;
data.problem.Dimension = 2;
% gname = 'random_poly_mesh_L1_n2_a0.9';
% gname = 'random_poly_mesh_L1_n4_a0.9';
% gname = 'random_poly_mesh_L1_n8_a0.9';    %ans paper
% gname = 'random_poly_mesh_L1_n16_a0.9';
% gname = 'random_poly_mesh_L1_n32_a0.9';
% gname = 'random_poly_mesh_L1_n64_a0.9';
% gname = 'shestakov_poly_mesh_L1_nc1_a0.1';
% gname = 'shestakov_poly_mesh_L1_nc2_a0.1';
% gname = 'shestakov_poly_mesh_L1_nc3_a0.1';
% gname = 'shestakov_poly_mesh_L1_nc4_a0.1';
% gname = 'shestakov_poly_mesh_L1_nc5_a0.1';      % this one
% gname = 'shestakov_poly_mesh_L1_nc6_a0.1';
% gname = 'shestakov_quad_L1_nc7_emb1_a0.2';
% gname = 'shestakov_quad_L1_nc7_emb2_a0.2';
% gname = 'shestakov_quad_L1_nc7_emb3_a0.2';
% gname = 'shestakov_quad_L1_nc7_emb4_a0.2';    %ans paper 
% gname = 'shestakov_quad_L1_nc7_emb5_a0.2';
% gname = 'shestakov_quad_L1_nc7_emb6_a0.2';
% gname = 'shestakov_quad_L1_nc7_emb7_a0.2';
% gname = 'smooth_quad_mesh_L1_nc7_emb1_a0.15';
% gname = 'smooth_quad_mesh_L1_nc7_emb2_a0.15';
% gname = 'smooth_quad_mesh_L1_nc7_emb3_a0.15';
% gname = 'smooth_quad_mesh_L1_nc7_emb4_a0.15';
% gname = 'smooth_quad_mesh_L1_nc7_emb5_a0.15';
% gname = 'smooth_quad_mesh_L1_nc7_emb6_a0.15';
% gname = 'smooth_quad_mesh_L1_nc7_emb7_a0.15';
% gname = 'smooth_poly_mesh_L1_n2_a0.15';
gname = 'smooth_poly_mesh_L1_n4_a0.15';
% gname = 'smooth_poly_mesh_L1_n8_a0.15';
% gname = 'smooth_poly_mesh_L1_n16_a0.15';
% gname = 'smooth_poly_mesh_L1_n32_a0.15'; % this one
% gname = 'z_mesh_quad_L1_n5_a0.05';
% gname = 'z_mesh_quad_L1_n6_a0.05';
% gname = 'z_mesh_quad_L1_n9_a0.05';
% gname = 'z_mesh_quad_L1_n20_a0.05';
% gname = 'z_mesh_quad_L1_n320_a0.05';
% gname = 'z_mesh_poly_L1_n9_a0.05';
% gname = 'z_mesh_poly_L1_n20_a0.05';
% gname = 'z_mesh_poly_L1_n40_a0.05';        % this one
% gname = 'z_mesh_poly_L1_n80_a0.05';
% gname = 'misha_quad_L1_n4';
load(strcat(glob.geom_path,gname,'.mat'));


% [x,y]=meshgrid(linspace(0,L,n+1),linspace(0,L,n+1));
% x=x(:);y=y(:);
% tri = delaunayTriangulation(x,y);
% geometry = GeneralGeometry(2, 'Delaunay', tri);

x=linspace(0,L,n+1);
y=linspace(0,L,n+1);
z=linspace(0,L,n+1);
% y=linspace(0,100,(n-1)/10+1);
% geometry = CartesianGeometry(1,x);
% geometry = CartesianGeometry(2,x,y);
% geometry = CartesianGeometry(3,x,y,z);

% geometry.extrude_mesh_2D_to_3D([0,.25*L,.5*L,.75*L,L]);

% geometry.set_face_flag_on_surface(2,[0,0;0,L]);
% geometry.set_face_flag_on_surface(3,[L,0;L,L]);
% geometry.set_face_flag_on_surface(2,[0,0;0,100]);
% geometry.set_face_flag_on_surface(3,[100,0;100,100]);
% geometry.set_face_flag_on_surface(2,[0,0,0;1,0,0;1,0,1;0,0,1]);
% geometry.set_face_flag_on_surface(3,[0,1,0;1,1,0;1,1,1;0,1,1]);

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
