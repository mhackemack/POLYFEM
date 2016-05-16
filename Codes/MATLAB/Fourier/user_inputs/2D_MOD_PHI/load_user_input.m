function data = load_user_input(n,L)
global glob
% outputs
data.Output.plotting_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry.type = 'cart';
data.geometry.x = L;
% log_xmin = 0; log_xmax = 0; xnum = 1;
% data.geometry.x = logspace(log_xmin, log_xmax, xnum);
data.geometry.dyz = [1];
data.geometry.ncellx = n;
data.geometry.ncelly = n;
data.geometry.ncellz = 1;
% mat regions
data.problem.NumberMaterials = 2;
% dy = 1/(n/2);
% for i=1:n/2
%     y0 = (i-1)*dy;
%     y1 = y0 + .5*dy;
%     y2 = y0 + dy;
%     data.geometry.mats(i).ID = 2;
%     data.geometry.mats(i).Region = [0,y1;1,y1;1,y2;0,y2];
% end
% data.geometry.mats(1).ID = 2;
% data.geometry.mats(1).Region = [0,0;.5,0;.5,.5;0,.5];
data.geometry.mats(1).ID = 2;
data.geometry.mats(1).Region = [0,.5;1,.5;1,1;0,1];
% fem
data.problem.refineMesh = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.FEMLumping = false;
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.TransportMethod = 'SI';
data.Neutronics.Transport.transportType = 'upwind';
% acceleration
data.Neutronics.PerformAcceleration = false;
data.Neutronics.DSAType = 'MIP';
data.Neutronics.AccelType = glob.Accel_WGS_DSA;
data.Neutronics.IP_Constant = 4;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [2,4,8,16];
data.Neutronics.Transport.PnOrder = 0;
% groups
data.Neutronics.numberEnergyGroups = 1;
% xs
% c = 0.9999; sigt = [1e-1;1e1];
% data.Neutronics.TotalXS = sigt;
% data.Neutronics.DiffusionXS = (1/3)./sigt;
% data.Neutronics.ScatteringXS = c*sigt;
% data.Neutronics.AbsorbXS = (1-c)*sigt;
data.Neutronics.ErrorShape = ones(data.problem.NumberMaterials,data.Neutronics.numberEnergyGroups);
% bcs
data.Neutronics.BCFlags = glob.Periodic;
data.Neutronics.BCVals = {0.0};
% average cross sections
% data.Neutronics.AveTotalXS      = data.Neutronics.TotalXS;
% data.Neutronics.AveDiffusionXS  = data.Neutronics.DiffusionXS;
% data.Neutronics.AveScatteringXS = data.Neutronics.ScatteringXS;
% data.Neutronics.AveAbsorbXS     = data.Neutronics.AbsorbXS;