function [data, geometry] = load_user_input()
global glob
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.problem.Path = 'Diffusion_MMS/Gauss_2D';
data.problem.Name = 'AMR_cart_PWLD_Irr=2_rtol=5.0';
data.problem.NumberMaterials = 1;
data.problem.problemType = 'SourceDriven';
data.problem.plotSolution = 0;
data.problem.saveSolution = 0;
data.problem.saveVTKSolution = 1;
% AMR Input Parameters
% ------------------------------------------------------------------------------
data.problem.refineMesh = 1;
data.problem.refinementLevels = 15;
data.problem.refinementTolerance = 0.5;
data.problem.AMRIrregularity = 2;
data.problem.projectSolution = 0;
data.problem.refinementType = 0; % 0 = err(c)/maxerr < c, 1 = numc/totalCells = c
% Neutronics Data
% ------------------------------------------------------------------------------
data.Neutronics.PowerLevel = 1.0; % only for eigenvalue problems
data.Neutronics.StartingSolution = 'zero';
data.Neutronics.transportMethod = 'Diffusion';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.SpatialMethod = 'MAXENT';
data.Neutronics.FEMDegree = 1;
data.Neutronics.numberEnergyGroups = 1;

%%%%%%%%%%%%%%%%%%%%%%
% Diffusion Properties
%%%%%%%%%%%%%%%%%%%%%%
data.Neutronics.Diffusion.MMS = true;
data.Neutronics.Diffusion.QuadOrder = 6;
data.Neutronics.Diffusion.ExtSource = cell(data.Neutronics.numberEnergyGroups,1);
data.Neutronics.Diffusion.ExactSolution = cell(data.Neutronics.numberEnergyGroups,1);
% Physical Properties
data.Neutronics.Diffusion.ScatteringXS = zeros(1,1,1);
% data.Neutronics.Diffusion.DiffXS = [1];
% data.Neutronics.Diffusion.TotalXS = [0];
% data.Neutronics.Diffusion.AbsorbXS = [0];
data.Neutronics.Diffusion.DiffXS = [1/3];
data.Neutronics.Diffusion.TotalXS = [1];
data.Neutronics.Diffusion.AbsorbXS = [1];
data.Neutronics.Diffusion.ScatteringXS(1,:,:) = [0.0];
data.Neutronics.Diffusion.FissionXS = [0.0];
data.Neutronics.Transport.NuBar = [0.0];
data.Neutronics.Diffusion.FissSpec = [0.0];
% data.Neutronics.Diffusion.ExtSource{1} = @rhs_func_quad_patch;
% data.Neutronics.Diffusion.ExactSolution{1} = @sol_func_quad_patch;
data.Neutronics.Diffusion.ExtSource{1} = @rhs_forcing_func_gauss2;
data.Neutronics.Diffusion.ExactSolution{1} = @mms_exact_solution_gauss2;
% Boundary Conditions
data.Neutronics.Diffusion.BCFlags = [glob.Function];
data.Neutronics.Diffusion.BCVals{1,1} = @mms_exact_solution_gauss2;

% Solver Input Parameters
% -----------------------
data.solver.absoluteTolerance = 1e-12;
data.solver.relativeTolerance = 1e-10;
data.solver.maxIterations = 2000;
data.solver.performNKA = 0;
data.solver.kyrlovSubspace = [];

% Geometry Data
% -------------
data.problem.Dimension = 2;
L = 1; ncells = 4;
% gname = 'random_poly_mesh_L1_n8_a0.9';    %ans paper
% gname = 'shestakov_poly_mesh_L1_nc5_a0.1';      % this one
% gname = 'shestakov_quad_L1_nc7_emb4_a0.2';    %ans paper
% gname = 'smooth_poly_mesh_L1_n8_a0.15';
% gname = 'z_mesh_poly_L1_n9_a0.05';
% gname = 'misha_quad_L1_n4';
% load(strcat(glob.geom_path,gname,'.mat'));

% xx=linspace(0,L,ncells+1);
% [x,y]=meshgrid(xx,xx);
% x=x(:);y=y(:);
% tri = delaunayTriangulation(x,y);
% geometry = GeneralGeometry(2, 'Delaunay', tri);
% geometry.extrude_mesh_2D_to_3D(xx)

x=linspace(0,L,ncells+1);
y=linspace(0,L,ncells+1);
z=linspace(0,L,ncells+1);
% geometry = CartesianGeometry(1,x);
geometry = CartesianGeometry(2,x,y);
% geometry = CartesianGeometry(3,x,y,z);

% geometry.set_face_flag_on_surface(2,[0,0;0,1]);
% geometry.set_face_flag_on_surface(3,[1,0;1,1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   MMS Function Listings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rhs_forcing_func_gauss2(xx)
x = xx(:,1);
y = xx(:,2);
Lx=1;
Ly=Lx;
x0=3*Lx/4;y0=x0;
varia=Lx^2/100;
c_diff = 1/3;
sigma_a = 1.0;

out = c_diff.*(1.0./Lx.^2.*1.0./Ly.^2.*x.*exp(-((x-x0).^2+(y-y0).^2)./varia).*...
      (Lx-x).*2.0e2+1.0./Lx.^2.*1.0./Ly.^2.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*...
      (Ly-y).*2.0e2-(1.0./Lx.^2.*1.0./Ly.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*...
      (x.*2.0-x0.*2.0).*(Ly-y).*2.0e2)./varia-(1.0./Lx.^2.*1.0./Ly.^2.*x.*y.*exp(-((x-x0).^2+...
      (y-y0).^2)./varia).*(y.*2.0-y0.*2.0).*(Lx-x).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*y.*...
      exp(-((x-x0).^2+(y-y0).^2)./varia).*(x.*2.0-x0.*2.0).*(Lx-x).*(Ly-y).*2.0e2)./varia+...
      (1.0./Lx.^2.*1.0./Ly.^2.*x.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(y.*2.0-y0.*2.0).*(Lx-x).*...
      (Ly-y).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*...
      (Lx-x).*(Ly-y).*4.0e2)./varia-1.0./Lx.^2.*1.0./Ly.^2.*1.0./varia.^2.*x.*y.*...
      exp(-((x-x0).^2+(y-y0).^2)./varia).*(x.*2.0-x0.*2.0).^2.*(Lx-x).*(Ly-y).*1.0e2-1.0./Lx.^2.*...
      1.0./Ly.^2.*1.0./varia.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(y.*2.0-y0.*2.0).^2.*(Lx-x).*...
      (Ly-y).*1.0e2)+1.0./Lx.^2.*1.0./Ly.^2.*sigma_a.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(Lx-x).*(Ly-y).*1.0e2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mms_exact_solution_gauss2(xx)
x = xx(:,1);
y = xx(:,2);
Lx=1;
Ly=Lx;
x0=3*Lx/4;y0=x0;
varia=Lx^2/100;

out = 1.0./Lx.^2.*1.0./Ly.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(Lx-x).*(Ly-y).*1.0e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rhs_forcing_func_quad2(xx)
x = xx(:,1);
y = xx(:,2);
Lx=1;
Ly=Lx;
c_diff = 1/3;
sigma_a = 1.0;

out = c_diff.*(x.*(Lx-x).*3.2e1+y.*(Ly-y).*3.2e1)+sigma_a.*x.*y.*(Lx-x).*(Ly-y).*1.6e1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mms_exact_solution_quad2(xx)
x = xx(:,1);
y = xx(:,2);
Lx=1;
Ly=Lx;

out = x.*y.*(Lx-x).*(Ly-y).*1.6e1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rhs_forcing_func_gauss3(xx)
x = xx(:,1);
y = xx(:,2);
z = xx(:,3);
Lx=1;
Ly=Lx;
Lz=Lx;
x0=3*Lx/4;y0=x0;z0=x0;
varia=Lx^2/10;
c_diff = 1/3;
sigma_a = 1.0;

out = c_diff.*(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(Lx-x).*(Ly-y).*2.0e2+1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(Lx-x).*(Lz-z).*2.0e2+1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(Ly-y).*(Lz-z).*2.0e2+(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(Lx-x).*(Ly-y).*(Lz-z).*6.0e2)./varia-(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(x.*2.0-x0.*2.0).*(Ly-y).*(Lz-z).*2.0e2)./varia-(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(y.*2.0-y0.*2.0).*(Lx-x).*(Lz-z).*2.0e2)./varia-(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(z.*2.0-z0.*2.0).*(Lx-x).*(Ly-y).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(x.*2.0-x0.*2.0).*(Lx-x).*(Ly-y).*(Lz-z).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(y.*2.0-y0.*2.0).*(Lx-x).*(Ly-y).*(Lz-z).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(z.*2.0-z0.*2.0).*(Lx-x).*(Ly-y).*(Lz-z).*2.0e2)./varia-1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*1.0./varia.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(x.*2.0-x0.*2.0).^2.*(Lx-x).*(Ly-y).*(Lz-z).*1.0e2-1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*1.0./varia.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(y.*2.0-y0.*2.0).^2.*(Lx-x).*(Ly-y).*(Lz-z).*1.0e2-1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*1.0./varia.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(z.*2.0-z0.*2.0).^2.*(Lx-x).*(Ly-y).*(Lz-z).*1.0e2)+1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*sigma_a.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(Lx-x).*(Ly-y).*(Lz-z).*1.0e2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mms_exact_solution_gauss3(xx)
x = xx(:,1);
y = xx(:,2);
z = xx(:,3);
Lx=1;
Ly=Lx;
Lz=Lx;
x0=3*Lx/4;y0=x0;z0=x0;
varia=Lx^2/10;

out = 1.0./Lx.^2.*1.0./Ly.^2.*1.0./Lz.^2.*x.*y.*z.*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)./varia).*(Lx-x).*(Ly-y).*(Lz-z).*1.0e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rhs_forcing_func_quad3(xx)
x = xx(:,1);
y = xx(:,2);
z = xx(:,3);
Lx=1;
Ly=Lx;
Lz=Lx;
c_diff = 1/3;
sigma_a = 1.0;

out = c_diff.*(x.*y.*(Lx-x).*(Ly-y).*2+x.*z.*(Lx-x).*(Lz-z).*2+y.*z.*(Ly-y).*(Lz-z).*2)+sigma_a.*x.*y.*z.*(Lx-x).*(Ly-y).*(Lz-z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mms_exact_solution_quad3(xx)
x = xx(:,1);
y = xx(:,2);
z = xx(:,3);
Lx=1;
Ly=Lx;
Lz=Lx;
out = x.*y.*z.*(Lx-x).*(Ly-y).*(Lz-z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rhs_func_quad_patch(~)
% x = xx(:,1); y = xx(:,2);
c_diff = 1;

out = c_diff.*-1.0e1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sol_func_quad_patch(xx)
x = xx(:,1); y = xx(:,2);

out = x.*-2.0+y.*6.0-x.*y.*3.0+x.^2+y.^2.*4.0+1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%