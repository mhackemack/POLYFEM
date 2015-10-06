function [data, geometry] = load_user_input()
global glob
% Problem Input Parameters
% ------------------------------------------------------------------------------
data.IO.Path = 'Transport/Homogeneous';
data.IO.Name = 'AMR_cart_Irr=3_rtol=0.2';
data.IO.PlotSolution = 0;
data.IO.SaveSolution = 0;
data.IO.SaveVTKSolution = 0;
data.IO.PrintIterationInfo = 1;
data.IO.PrintConstructionInfo = 1;
% AMR Input Parameters
% ------------------------------------------------------------------------------
data.AMR.RefineMesh = 0;
% MMS Input Parameters
% ------------------------------------------------------------------------------
data.MMS.PerformMMS = 0;
% Overall Problem Data
% ------------------------------------------------------------------------------
data.problem.NumberMaterials = 1;
data.problem.ProblemType = 'SourceDriven';
data.problem.PowerLevel = 1.0;
data.problem.TransportMethod = 'Transport';
data.problem.FEMType = 'DFEM';
data.problem.SpatialMethod = 'PWLD';
data.problem.FEMDegree = 1;
% Transport Properties
% ------------------------------------------------------------------------------
% Tranpsort Type Properties - most of this only applies to hybrid transport
data.Transport.PerformAcceleration = true;
data.Transport.PnOrder = 0;
data.Transport.XSID = 1; data.Transport.QuadID = 1;
data.Transport.TransportType = 'upwind';
data.Transport.PerformSweeps = 0;
data.Transport.VisualizeSweeping = 0;
data.Transport.StabilizationMethod = 'EGDG';
data.Transport.FluxStabilization = 2.0;
data.Transport.CurrentStabilization = 1.0;
% Quadrature Properties
data.Quadrature(1).PnOrder = data.Transport.PnOrder;
data.Quadrature(1).AngleAggregation = 'all';
data.Quadrature(1).QuadType = 'LS';
data.Quadrature(1).SnLevels = 4;
data.Quadrature(1).PolarLevels = 4;
data.Quadrature(1).AzimuthalLevels = 4;
% Flux Properties
data.Fluxes.StartingSolution = 'zero';
% Retrieve All Physical Properties
data = get_69G_Graphite_XS(data, data.Transport.PnOrder);
data.XS(1).ExtSource = zeros(1,data.Groups.NumberEnergyGroups);
data.XS(1).ExtSource(1,1) = 1.0;
data.XS(1).BCFlags = [glob.Vacuum]; data.XS(1).BCVals{1} = 0;
data.XS(2).BCFlags = [glob.Vacuum]; data.XS(2).BCVals{1} = 0;
% Construct Group Set Information
data.Groups.NumberGroupSets = data.Groups.NumberEnergyGroups;
data.Groups.GroupSets = cell(data.Groups.NumberGroupSets,1);
data.Groups.GroupSetUpscattering = [false(length(data.Groups.FastGroups),1);true(length(data.Groups.ThermalGroups),1)];
for g=1:data.Groups.NumberGroupSets
    data.Groups.GroupSets{g} = g;
end
data = collapse_two_grid_xs(data);
% Acceleration Properties
% ------------------------------------------------------------------------------
data.Acceleration.WGSAccelerationBool = false(data.Groups.NumberGroupSets,1);
data.Acceleration.AGSAccelerationBool = true;
data.Acceleration.WGSAccelerationResidual = false(data.Groups.NumberGroupSets,1);
data.Acceleration.AGSAccelerationResidual = false;
data.Acceleration.WGSAccelerationID = zeros(data.Groups.NumberGroupSets,1);
data.Acceleration.AGSAccelerationID = 1;
data.Acceleration.Info(1).AccelerationType = glob.Accel_AGS_TG;
data.Acceleration.Info(1).DiscretizationType = glob.Accel_DSA_MIP;
data.Acceleration.Info(1).IP_Constant = 4;
data.Acceleration.Info(1).Groups = data.Groups.ThermalGroups;
data.Acceleration.Info(1).Moments = 1;
data.Acceleration.Info(1).XSID = 2;
% Solver Input Parameters
% ------------------------------------------------------------------------------
data.solver.AGSMaxIterations = 500000;
data.solver.WGSMaxIterations = 1e2*ones(data.Groups.NumberGroupSets,1);
data.solver.AGSRelativeTolerance = 1e-4;
data.solver.WGSRelativeTolerance = 1e-6*ones(data.Groups.NumberGroupSets,1);
data.solver.AGSAbsoluteTolerance = 1e-4;
data.solver.WGSAbsoluteTolerance = 1e-6*ones(data.Groups.NumberGroupSets,1);
% Geometry Data
% ------------------------------------------------------------------------------
data.problem.Dimension = 2;
L = 1e1; ncells = 10;
% L = 30*6.953164954422388e-02; ncells = 30;

x=linspace(0,L,ncells+1);
y=linspace(0,L,ncells+1);
z=linspace(0,L,ncells+1);
% geometry = CartesianGeometry(1,x);
geometry = CartesianGeometry(2,x,y);
% geometry = CartesianGeometry(3,x,y,z);

% geometry.set_face_flag_on_surface(2,0.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = collapse_two_grid_xs(data)
nm = data.problem.NumberMaterials;
grps = data.Groups.ThermalGroups; ngrps = length(grps);
txs = data.XS(1).TotalXS(:,grps);
sxs = data.XS(1).ScatteringXS(:,grps,grps,1);
data.Acceleration.Info(1).ErrorShape = zeros(nm,ngrps);
% Generate Eigenshape for energy collaps
y = cell(nm,1);
for m=1:nm
    A = (diag(txs(m,:)) - tril(squeeze(sxs(m,:,:))))\triu(squeeze(sxs(m,:,:)),1);
    [y{m},~,~] = power_method(A,ones(ngrps,1),2000,1e-15);
    y{m} = y{m} / sum(y{m});
    data.Acceleration.Info(1).ErrorShape(m,:) = y{m};
end
% Populate Acceleration XS
for m=1:nm
    data.XS(2).TotalXS(m) = 0;
    data.XS(2).DiffXS(m) = 0;
    data.XS(2).AbsorbXS(m) = 0;
    for g=1:ngrps
        data.XS(2).TotalXS(m) = data.XS(2).TotalXS(m) + y{m}(g)*txs(m,g);
        data.XS(2).DiffXS(m) = data.XS(2).DiffXS(m) + y{m}(g)/txs(m,g)/3;
        data.XS(2).AbsorbXS(m) = data.XS(2).AbsorbXS(m) + y{m}(g)*txs(m,g);
        for gg=1:ngrps
            data.XS(2).AbsorbXS(m) = data.XS(2).AbsorbXS(m) - sxs(m,g,gg,1)*y{m}(gg);
        end
    end
end
data.XS(2).ScatteringXS = data.XS(1).ScatteringXS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%