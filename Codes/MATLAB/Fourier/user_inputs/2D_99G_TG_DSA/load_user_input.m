function data = load_user_input()
global glob
% outputs
data.Output.plotting_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry.type = 'cart';
log_xmin = -1; log_xmax = -1; xnum = 1;
data.geometry.x = logspace(log_xmin, log_xmax, xnum);
data.geometry.dyz = [1];
data.geometry.ncellx = 2;
data.geometry.ncelly = 2;
data.geometry.ncellz = 1;
% mat regions
data.problem.NumberMaterials = 2;
data.geometry.mats = [];
data.geometry.mats(1).ID = 2;
data.geometry.mats(1).Region = [0,0;.5,0;.5,.5;0,.5];
data.geometry.mats(2).ID = 2;
data.geometry.mats(2).Region = [.5,.5;1,.5;1,1;.5,1];
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
data.Neutronics.AccelType = glob.Accel_AGS_TG;
data.Neutronics.IP_Constant = 4;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [4];
data.Neutronics.Transport.PnOrder = 0;
data.Transport.PnOrder = 0;
% groups
data.Neutronics.numberEnergyGroups = 57;
data.Neutronics.ThermalGroups = 43:99;
% bcs
data.Neutronics.BCFlags = glob.Periodic;
data.Neutronics.BCVals = {0.0};
% Build cross sections
% ------------------------------------------------------------------------------
% Rev1 changes
data.Groups.NumberEnergyGroups = 99;
data.Acceleration.Info.Groups = 1:data.Neutronics.numberEnergyGroups;
data.Acceleration.Info.AccelerationType = data.Neutronics.AccelType;
% graphite
% data = add_xs_component_contribution(data, 1, 2, 'graphite_99G', 8.5238E-2);
% HDPE
data = add_xs_component_contribution(data, 1, 2, 'PolyH1_99G', 8.1570E-2);
data = add_xs_component_contribution(data, 1, 2, 'FG_CNat_99G', 4.0787E-2);
% BHPDE
% data = add_xs_component_contribution(data, 1, 3, 'PolyH1_99G', 5.0859E-2);
% data = add_xs_component_contribution(data, 1, 3, 'FG_CNat_99G', 2.5429E-2);
% data = add_xs_component_contribution(data, 1, 3, 'B10_99G', 6.6256E-3);
% data = add_xs_component_contribution(data, 1, 3, 'B11_99G', 2.6669E-2);
% Air
data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 7.4906E-9);
data = add_xs_component_contribution(data, 1, 1, 'N14_99G', 3.9123E-5);
data = add_xs_component_contribution(data, 1, 1, 'O16_99G', 1.0511E-5);
data = add_xs_component_contribution(data, 1, 1, 'Ar40_99G', 2.3297E-7);
% Wood
% data = add_xs_component_contribution(data, 1, 1, 'FG_H1_99G', 2.0752E-2);
% data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 1.4520E-2);
% data = add_xs_component_contribution(data, 1, 1, 'O16_99G', 1.0376E-2);
% Restrict Energy Groups
data.XS.TotalXS = data.XS.TotalXS(:,data.Neutronics.ThermalGroups);
data.XS.AbsorbXS = data.XS.AbsorbXS(:,data.Neutronics.ThermalGroups);
data.XS.ScatteringXS = data.XS.ScatteringXS(:,data.Neutronics.ThermalGroups,data.Neutronics.ThermalGroups,:);
data.Neutronics.DiffusionXS = [];
data.Neutronics.TotalXS = data.XS.TotalXS;
data.Neutronics.AbsorbXS = data.XS.AbsorbXS;
data.Neutronics.ScatteringXS = data.XS.ScatteringXS;
% Collapse energy groups
data = collapse_tg_xs(data,1,2,1);
% Average cross sections
data.Neutronics.AveTotalXS = [];
data.Neutronics.AveScatteringXS = [];
data.Neutronics.AveDiffusionXS = data.XS(2).DiffXS;
data.Neutronics.AveAbsorbXS = data.XS(2).AbsorbXS;
data.Neutronics.ErrorShape = data.Acceleration.Info.ErrorShape;
data = rmfield(data,'XS');
data = get_analytical_solutions(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_analytical_solutions(data)
nm = data.problem.NumberMaterials;
ng = data.Neutronics.numberEnergyGroups;
for m=1:nm
    T = diag(data.Neutronics.TotalXS(m,:));
    S = squeeze(data.Neutronics.ScatteringXS(m,:,:));
    Sd = tril(S,0); Su = triu(S,1);
    F = (T-Sd)\Su; I = eye(ng); FF = Su*(F - I);
    [V,D] = eig(F); D=(diag(D));
    data.AnalyticalUnacceleratedEigenvalue(m,:) = D;
    data.AnalyticalUnacceleratedEigenvector{m} = V;
    [data.AnalyticalMaxUnacceleratedEigenvalue(m),ind] = max(abs(D));
    data.AnalyticalMaxUnacceleratedEigenvector(m,:) = V(:,ind);
    V = data.Neutronics.ErrorShape(m,:)';
    L = F + V*(sum((T-Sd-Su)*V))^(-1)*sum(FF,1);
    [V,D] = eig(L); D=(diag(D));
    data.AnalyticalAcceleratedEigenvalue(m,:) = D;
    data.AnalyticalAcceleratedEigenvector{m} = V;
    [data.AnalyticalMaxAcceleratedEigenvalue(m),ind] = max(abs(D));
    data.AnalyticalMaxAcceleratedEigenvector(m,:) = V(:,ind);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%