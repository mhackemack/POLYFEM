function data = load_user_input()
global glob
% outputs
data.Output.plotting_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 1;
data.geometry.type = 'cart';
data.geometry.x = [1e0];
% data.geometry.x = [1e0,1e-1,1e-2,1e-3];
% log_xmin = -4; log_xmax = 0; xnum = 5;
% data.geometry.x = logspace(log_xmin, log_xmax, xnum);
data.geometry.dyz = [1];
data.geometry.ncellx = 1;
data.geometry.ncelly = 1;
data.geometry.ncellz = 1;
% mat regions
data.problem.NumberMaterials = 1;
% data.geometry.mats(1).ID = 1;
% data.geometry.mats(1).Region = [0.5;1];
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
data.Neutronics.PerformAcceleration = 0;
data.Neutronics.DSAType = 'MIP';
data.Neutronics.AccelType = glob.Accel_WGS_DSA;
data.Neutronics.IP_Constant = 4;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [8];
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
% data = add_xs_component_contribution(data, 1, 1, 'graphite_99G', 8.5238E-2);
% HDPE
% data = add_xs_component_contribution(data, 1, 1, 'PolyH1_99G', 8.1570E-2);
% data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 4.0787E-2);
% BHPDE
% data = add_xs_component_contribution(data, 1, 1, 'PolyH1_99G', 5.0859E-2);
% data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 2.5429E-2);
% data = add_xs_component_contribution(data, 1, 1, 'B10_99G', 6.6256E-3);
% data = add_xs_component_contribution(data, 1, 1, 'B11_99G', 2.6669E-2);
% Air
% data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 7.4906E-9);
% data = add_xs_component_contribution(data, 1, 1, 'N14_99G', 3.9123E-5);
% data = add_xs_component_contribution(data, 1, 1, 'O16_99G', 1.0511E-5);
% data = add_xs_component_contribution(data, 1, 1, 'Ar40_99G', 2.3297E-7);
% Wood
% data = add_xs_component_contribution(data, 1, 1, 'FG_H1_99G', 2.0752E-2);
% data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 1.4520E-2);
% data = add_xs_component_contribution(data, 1, 1, 'O16_99G', 1.0376E-2);
% AmBe
% data = add_xs_component_contribution(data, 1, 1, 'Am241_99G', 1.1649E-3);
% data = add_xs_component_contribution(data, 1, 1, 'Be9_99G', 1.9077E-1);
% data = add_xs_component_contribution(data, 1, 1, 'O16_99G', 2.3298E-3);
% Steel
data = add_xs_component_contribution(data, 1, 1, 'Cr52_99G', 1.7428E-2);
data = add_xs_component_contribution(data, 1, 1, 'Mn55_99G', 1.7363E-3);
data = add_xs_component_contribution(data, 1, 1, 'Fe56_99G', 5.9358E-2);
data = add_xs_component_contribution(data, 1, 1, 'Ni58_99G', 7.7199E-3);
% Boral
% data = add_xs_component_contribution(data, 1, 1, 'Al_99G', 3.8193E-2);
% data = add_xs_component_contribution(data, 1, 1, 'B10_99G', 7.1036E-3);
% data = add_xs_component_contribution(data, 1, 1, 'B11_99G', 2.8593E-2);
% data = add_xs_component_contribution(data, 1, 1, 'FG_CNat_99G', 8.9241E-3);
% BF3
% data = add_xs_component_contribution(data, 1, 1, 'B10_99G', 6.4458E-6);
% data = add_xs_component_contribution(data, 1, 1, 'B11_99G', 2.6858E-7);
% data = add_xs_component_contribution(data, 1, 1, 'F19_99G', 2.0143E-5);
% Restrict Energy Groups
data.XS.TotalXS = data.XS.TotalXS(:,data.Neutronics.ThermalGroups);
data.XS.AbsorbXS = data.XS.AbsorbXS(:,data.Neutronics.ThermalGroups);
data.XS.ScatteringXS = data.XS.ScatteringXS(:,data.Neutronics.ThermalGroups,data.Neutronics.ThermalGroups,:);
data.Neutronics.DiffusionXS = [];
data.Neutronics.TotalXS = data.XS.TotalXS;
data.Neutronics.AbsorbXS = data.XS.AbsorbXS;
data.Neutronics.ScatteringXS = data.XS.ScatteringXS;
% Collapse energy groups
data = collapse_jacobi_xs(data,1,2,1);
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
    F = T\S; I = eye(ng); FF = S*(F - I);
    [V,D] = eig(F); D=(diag(D));
    data.AnalyticalUnacceleratedEigenvalue(m,:) = D;
    data.AnalyticalUnacceleratedEigenvector{m} = V;
    [data.AnalyticalMaxUnacceleratedEigenvalue(m),ind] = max(abs(D));
    data.AnalyticalMaxUnacceleratedEigenvector(m,:) = V(:,ind);
    V = data.Neutronics.ErrorShape(m,:)';
    L = F + V*(sum((T-S)*V))^(-1)*sum(FF,1);
    [V,D] = eig(L); D=(diag(D));
    data.AnalyticalAcceleratedEigenvalue(m,:) = D;
    data.AnalyticalAcceleratedEigenvector{m} = V;
    [data.AnalyticalMaxAcceleratedEigenvalue(m),ind] = max(abs(D));
    data.AnalyticalMaxAcceleratedEigenvector(m,:) = V(:,ind);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%