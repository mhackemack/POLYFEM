function data = load_user_input()
global glob
% outputs
data.Output.plotting_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry.type = 'cart';
log_xmin = 0; log_xmax = 0; xnum = 1;
data.geometry.x = logspace(log_xmin, log_xmax, xnum);
data.geometry.dyz = [1];
data.geometry.ncellx = 1;
data.geometry.ncelly = 1;
data.geometry.ncellz = 1;
% mat regions
data.problem.NumberMaterials = 1;
data.geometry.mats = [];
% data.geometry.mats(1).ID = 2;
% data.geometry.mats(1).Region = [0,0;.5,0;.5,.5;0,.5];
% data.geometry.mats(2).ID = 2;
% data.geometry.mats(2).Region = [.5,.5;1,.5;1,1;.5,1];
% fem
data.problem.refineMesh = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.FEMLumping = false;
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.TransportMethod = 'SI';
data.Neutronics.Transport.transportType = 'upwind';
% acceleration
data.Neutronics.DSAType = 'MIP';
data.Neutronics.AccelType = glob.Accel_WGS_DSA;
data.Neutronics.IP_Constant = 4;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [4];
data.Neutronics.Transport.PnOrder = 0;
% groups
data.Neutronics.numberEnergyGroups = 57;
data.Neutronics.ThermalGroups = 43:99;
% bcs
data.Neutronics.BCFlags = glob.Periodic;
data.Neutronics.BCVals = 0.0;
% Build graphite xs
% ------------------------------------------------------------------------------
xs_dir = [glob.xs_dir,'graphite_99G/'];
% Allocate Memory
ng = length(data.Neutronics.ThermalGroups);
data.Neutronics.TotalXS = zeros(1,ng);
data.Neutronics.ScatteringXS = zeros(1,ng,ng);
% Total XS
M = open([xs_dir,'MT_1.mat']);
data.Neutronics.TotalXS = M.mat(data.Neutronics.ThermalGroups)';
data.Neutronics.DiffusionXS = (1/3)./data.Neutronics.TotalXS;
% Scattering XS
M = open([xs_dir,'MT_2500.mat']);
data.Neutronics.ScatteringXS(1,:,:) = M.mat(data.Neutronics.ThermalGroups,data.Neutronics.ThermalGroups,1);
% Remaining XS
data.Neutronics.AbsorbXS = [];
% Collapse two-grid spectrum
T = diag(data.Neutronics.TotalXS);
S = squeeze(data.Neutronics.ScatteringXS(1,:,:));
A = (T - tril(S))\triu(S,1);
[y,~,~] = power_method(A,ones(ng,1),2000,1e-15);
data.Neutronics.EnergyShape = (y / sum(y))';
% Average cross sections
data.Neutronics.AveTotalXS = data.Neutronics.EnergyShape*data.Neutronics.TotalXS';
data.Neutronics.AveDiffusionXS = data.Neutronics.EnergyShape*data.Neutronics.DiffusionXS';
data.Neutronics.AveAbsorbXS = data.Neutronics.EnergyShape*data.Neutronics.TotalXS' - sum(S*y);
data.Neutronics.AveScatteringXS = [];