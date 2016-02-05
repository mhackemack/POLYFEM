function data = load_user_input()
global glob
% outputs
data.Output.plotting_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry.type = 'tri';
log_xmin = 0; log_xmax = 0; xnum = 1;
data.geometry.x = logspace(log_xmin, log_xmax, xnum);
data.geometry.dyz = [1];
data.geometry.ncellx = 1;
data.geometry.ncelly = 1;
data.geometry.ncellz = 1;
% mat regions
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
data.Neutronics.DSAType = 'MIP';
data.Neutronics.AccelType = glob.Accel_DSA_MIP;
data.Neutronics.IP_Constant = 4;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [4];
data.Neutronics.Transport.PnOrder = 0;
% groups
data.Neutronics.numberEnergyGroups = 1;
% xs
c = 0.9999; sigt = 1.0;
data.Neutronics.Transport.TotalXS = sigt;
data.Neutronics.Diffusion.DiffusionXS = (1/3)./sigt;
data.Neutronics.Transport.ScatteringXS = c*sigt;
data.Neutronics.Diffusion.AbsorbXS = (1-c)*sigt;
% bcs
data.Neutronics.Transport.BCFlags = glob.Periodic;
data.Neutronics.Transport.BCVals = 0.0;
% average cross sections
data.Neutronics.Transport.AveTotalXS = data.Neutronics.Transport.TotalXS;
data.Neutronics.Diffusion.AveDiffusionXS = data.Neutronics.Diffusion.DiffusionXS;
data.Neutronics.Transport.AveScatteringXS = data.Neutronics.Transport.ScatteringXS;
data.Neutronics.Diffusion.AveAbsorbXS = data.Neutronics.Diffusion.AbsorbXS;