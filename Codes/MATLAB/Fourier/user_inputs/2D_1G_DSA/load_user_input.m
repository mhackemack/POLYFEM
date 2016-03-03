function data = load_user_input()
global glob
% outputs
data.Output.plotting_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry.type = 'cart';
% data.geometry.x = [1e-1,1e0,1e1,1e2,1e3];
data.geometry.x = [1e0];
% log_xmin = 0; log_xmax = 0; xnum = 1;
% log_xmin = -3; log_xmax = 3; xnum = 241;
% data.geometry.x = logspace(log_xmin, log_xmax, xnum);
data.geometry.dyz = [1];
data.geometry.ncellx = 1;
data.geometry.ncelly = 1;
data.geometry.ncellz = 1;
% mat regions
data.problem.NumberMaterials = 1;
data.geometry.mats = [];
% fem
data.problem.refineMesh = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.FEMLumping = false;
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.TransportMethod = 'SI';
data.Neutronics.Transport.transportType = 'upwind';
% acceleration
data.Neutronics.PerformAcceleration = 1;
data.Neutronics.DSAType = 'MIP';
data.Neutronics.AccelType = glob.Accel_WGS_DSA;
data.Neutronics.IP_Constant = 4;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [4];
data.Neutronics.Transport.PnOrder = 0;
% groups
data.Neutronics.numberEnergyGroups = 1;
% xs
c = 0.9999; sigt = 1e-2;
data.Neutronics.TotalXS = sigt;
data.Neutronics.DiffusionXS = (1/3)./sigt;
data.Neutronics.ScatteringXS = c*sigt;
data.Neutronics.AbsorbXS = (1-c)*sigt;
data.Neutronics.ErrorShape = ones(data.problem.NumberMaterials,data.Neutronics.numberEnergyGroups);
% bcs
data.Neutronics.BCFlags = glob.Periodic;
data.Neutronics.BCVals = {0.0};
% average cross sections
data.Neutronics.AveTotalXS      = data.Neutronics.TotalXS;
data.Neutronics.AveDiffusionXS  = data.Neutronics.DiffusionXS;
data.Neutronics.AveScatteringXS = data.Neutronics.ScatteringXS;
data.Neutronics.AveAbsorbXS     = data.Neutronics.AbsorbXS;