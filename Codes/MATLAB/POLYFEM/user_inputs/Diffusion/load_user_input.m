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
data.problem.projectSolution = 1;
data.problem.plotSolution = 0;
data.problem.saveSolution = 0;

% Neutronics Data
% ---------------
data.Neutronics.PowerLevel = 1.0;
data.Neutronics.StartingSolution = 'random';
data.Neutronics.transportMethod = 'Diffusion';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMDegree = 1;
data.Neutronics.numberEnergyGroups = 1;

% Diffusion Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.Neutroncis.IP_Constant = 4;
data.Neutronics.Diffusion.MMS = false;
% Physical Properties
data.Neutronics.Diffusion.ScatteringXS = zeros(1,1,1);
% data.Neutronics.Diffusion.DiffXS = [1/3];
% data.Neutronics.Diffusion.TotalXS = [1];
data.Neutronics.Diffusion.DiffXS = [2];
data.Neutronics.Diffusion.TotalXS = [0];
data.Neutronics.Diffusion.AbsorbXS = [0];
data.Neutronics.Diffusion.ScatteringXS(1,:,:) = [0.0];
data.Neutronics.Diffusion.FissionXS = [0.0];
data.Neutronics.Transport.NuBar = [0.0];
data.Neutronics.Diffusion.FissSpec = [0.0];
data.Neutronics.Diffusion.ExtSource = [0.0];
% Boundary Conditions
% data.Neutronics.Diffusion.BCFlags = [glob.Dirichlet];
% data.Neutronics.Diffusion.BCVals = [0.0];
data.Neutronics.Diffusion.BCFlags = [glob.Neumann; glob.Robin; glob.Robin];
data.Neutronics.Diffusion.BCVals = [0.0;0.0;9.0];

% Solver Input Parameters
% -----------------------
data.solver.absoluteTolerance = 1e-12;
data.solver.relativeTolerance = 1e-10;
data.solver.maxIterations = 2000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Geometry Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% gname = 'smooth_poly_mesh_L1_n4_a0.15';
gname = 'smooth_poly_mesh_L1_n8_a0.15';
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

n=5;
% [x,y]=meshgrid(linspace(0,1,n),linspace(0,1,n));
% x=x(:);y=y(:);
% tri = delaunayTriangulation(x,y);
% geometry = GeneralGeometry(2, 'Delaunay', tri);

x=linspace(0,1,n);
y=linspace(0,1,n);
z=linspace(0,1,n);
% y=linspace(0,100,(n-1)/10+1);
% geometry = CartesianGeometry(1,x);
% geometry = CartesianGeometry(2,x,y);
% geometry = CartesianGeometry(3,x,y,z);

% geometry.extrude_mesh_2D_to_3D([0,.25,.5,.75,1]);

geometry.set_face_flag_on_surface(2,[0,0;0,1]);
geometry.set_face_flag_on_surface(3,[1,0;1,1]);
% geometry.set_face_flag_on_surface(2,[0,0;0,100]);
% geometry.set_face_flag_on_surface(3,[100,0;100,100]);
% geometry.set_face_flag_on_surface(2,[0,0,0;1,0,0;1,0,1;0,0,1]);
% geometry.set_face_flag_on_surface(3,[0,1,0;1,1,0;1,1,1;0,1,1]);



