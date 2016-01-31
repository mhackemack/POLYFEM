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
data.problem.ProblemType = 'Eigenvalue';
data.problem.PowerLevel = 1.0;
data.problem.TransportMethod = 'Transport';
data.problem.FEMType = 'DFEM';
data.problem.SpatialMethod = 'PWLD';
data.problem.FEMLumping = false;
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
% Flux Properties
data.Fluxes.StartingSolution = 'one';