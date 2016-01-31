%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execute Problem Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, geometry, DoF, FE] = execute_problem_Rev1(data, geometry)
global glob
% Generate problem execution schedule
% ------------------------------------------------------------------------------
pcall = get_execution_function_handle(data);
geometry = createAMRMesh(data, geometry);
[mms, deg] = get_mms_information(data);
% ------------------------------------------------------------------------------
% Perform Initial Solution Calculation (terminates if no mesh refinement)
% ------------------------------------------------------------------------------
% Get Degree of Freedom Handler - gets remade after a mesh refinement
DoF = DoFHandler(geometry, data.problem.FEMDegree, data.problem.FEMType, data.problem.DoFType);
% Get FEHandler based on MMS - this is necessary so that FEHandler
% generates Quadrature Sets for the Polygonal/Polyhedral cells
FE = FEHandler(geometry, DoF, data.problem.SpatialMethod, data.problem.FEMLumping, data.problem.FEMVolumeBools, data.problem.FEMSurfaceBools, mms, deg);
% Allocate Solution Space - gets reallocated after a mesh refinement
data = prepare_problem_execution(data, geometry);
% data = solution_allocation(data, geometry, DoF);
% Solve the neutronics problem (either source-driven or k-eigenvalue)
data = pcall(data, geometry, DoF, FE);
% sol.CellVertexNumbers = geometry.CellVertexNumbers;
% Calculate MMS Error if necessary
% if mms
%     sol.MMS_error = calculate_MMS_error(data, geometry, DoF, FE, sol.flux);
% end
% Plot solution for viewer to see
% if data.IO.PlotSolution
%     glob.plot_counter = glob.plot_counter + 1;
%     figure(glob.plot_counter)
%     plot_solution(geometry,DoF,FE,sol.flux)
% end
% Save off output objects if necessary
% if data.IO.SaveSolution
%     f_name = ['outputs/',data.IO.Path,'/',data.IO.Name];
%     o_str = '_0';
%     save([f_name,'_data',o_str,'.mat'], 'data');
%     save([f_name,'_geometry',o_str,'.mat'], 'geometry');
%     save([f_name,'_DoF',o_str,'.mat'], 'DoF');
%     save([f_name,'_FE',o_str,'.mat'], 'FE');
%     save([f_name,'_sol',o_str,'.mat'], 'sol');
% end
% Save off solution to VTK output for viewing in Visit/Paraview
% if data.IO.SaveVTKSolution
%     f_name = ['outputs/',data.IO.Path,'/',data.IO.Name];
%     sol_flux = cell(data.Groups.NumberEnergyGroups,1);
%     sol_name = cell(data.Groups.NumberEnergyGroups,1);
%     for i=1:data.Groups.NumberEnergyGroups
%         sol_flux{i} = sol.flux{i,1};
%         sol_name{i} = ['flux_g',num2str(i)];
%     end
%     write_output_to_vtk_rev2(f_name, data, geometry, DoF, sol_flux, sol_name);
% end
% ------------------------------------------------------------------------------
% Perform Mesh Refinement Calculations
% ------------------------------------------------------------------------------
% if data.AMR.RefineMesh && data.AMR.RefinementLevels > 0
%     % Retrieve some quick info first
%     tsol = sol; clear sol;
%     sol{1} = tsol;
%     % Loop through refinement levels
%     for iter=1:data.AMR.RefinementLevels
%         o_str = ['_',num2str(iter)];
%         r = iter + 1;
%         % Save old objects for flux interpolation if necessary
%         if data.AMR.ProjectSolution
%             data0 = data;
%             DoF0  = DoF;
%         end
%         % Form new mesh/DoFHandler/FEHandler
%         refine_problem_mesh(data, geometry, DoF, FE, sol{r-1}.flux); % This works because of pass-by-reference
%         DoF = DoFHandler(geometry, data.problem.FEMDegree, data.problem.FEMType, data.problem.DoFType);
%         FE = FEHandler(geometry, DoF, data.problem.SpatialMethod, data.Neutronics.FEMLumping, data.problem.FEMVolumeBools, data.problem.FEMSurfaceBools, mms, deg);
%         data = prepare_problem_execution(data, geometry);
%         [data, sol{r}] = solution_allocation(data, geometry, DoF);
%         % Interpolate flux solutions to next refinement level
%         if data.problem.projectSolution
%             [data,sol{r}.flux] = interpolate_ref_flux(data0,data,geometry,DoF0,DoF,sol{r-1}.flux,sol{r}.flux);
%         end
%         % Compute new flux solution
%         [data, tsol] = pcall(data, geometry, DoF, FE, sol{r});
%         sol{r} = tsol;
%         sol{r}.CellVertexNumbers = geometry.CellVertexNumbers;
%         % Calculate MMS Error if necessary
%         if mms
%             sol{r}.MMS_error = calculate_MMS_error(data, geometry, DoF, FE, sol{r}.flux);
%         end
%         % Save off output objects if necessary
%         if data.problem.SaveSolution
%             save([f_name,'_data',o_str,'.mat'], 'data');
%             save([f_name,'_geometry',o_str,'.mat'], 'geometry');
%             save([f_name,'_DoF',o_str,'.mat'], 'DoF');
%             save([f_name,'_FE',o_str,'.mat'], 'FE');
%             % Solution Vector
%             tsol = sol; sol = tsol{end};
%             save([f_name,'_sol',o_str,'.mat'], 'sol');
%             sol = tsol;
%         end
%         % Save off solution to VTK output for viewing in Visit/Paraview
%         if data.problem.SaveVTKSolution
%             f_name = ['outputs/',data.problem.Path,'/',data.problem.Name];
%             sol_flux = cell(data.Neutronics.numberEnergyGroups,1);
%             sol_name = cell(data.Neutronics.numberEnergyGroups,1);
%             for i=1:data.Neutronics.numberEnergyGroups
%                 sol_flux{i} = sol{r}.flux{i,1};
%                 sol_name{i} = ['flux_g',num2str(i)];
%             end
%             write_output_to_vtk_rev2(f_name, data, geometry, DoF, sol_flux, sol_name);
%         end
%     end
% end
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Procedure for source-driven transport problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = perform_transport_source_driven_problem(data, mesh, DoF, FE)
% Get Necessary Function Handles
% ------------------------------------------------------------------------------
pcall = @solve_linear_iteration_transport;
% Get some additional information
% ------------------------------------------------------------------------------
xsid = data.Transport.XSID;
qid  = data.Transport.QuadID;
% Set fission and external source contributions
% ------------------------------------------------------------------------------
ng = data.Groups.NumberEnergyGroups;
data.Sources.FissionSource = SetFissionSource_Transport(data.XS(xsid),1:ng,1:ng,data.Fluxes.Phi,mesh,DoF,FE);
data.Sources.ExtSource = SetExtSource_Transport(data,xsid,qid,1:ng,mesh,DoF,FE);
% Perform single source-driven linear iteration
% ------------------------------------------------------------------------------
data = pcall(data,xsid,qid,mesh,DoF,FE);
% Process some iteration/output information
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Procedure for source-driven diffusion problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = perform_diffusion_source_driven_problem(data, mesh, DoF, FE)
% Get Linear Iteration Function Handle
% ------------------------------------------------------------------------------
pcall = @solve_linear_iteration_diffusion;
% Get some additional information
% ------------------------------------------------------------------------------
trans_xsid = data.Diffusion.XSID;
% Perform single source-driven linear iteration
% ------------------------------------------------------------------------------
data = pcall(data,trans_xsid,mesh,DoF,FE);
% Process some iteration/output information
% ------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Iterative Procedure for k-eigenvalue transport problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = perform_transport_eignevalue_problem(data, mesh, DoF, FE)
% Get some transport information
% ------------------------------------------------------------------------------
xsid = data.Transport.XSID;
qid  = data.Transport.QuadID;
ng = data.Groups.NumberEnergyGroups;
% Set initial source terms
% ------------------------------------------------------------------------------
KeffOldPhi = data.Fluxes.Phi;
data.Sources.ExtSource = SetExtSource_Transport(data,xsid,qid,1:ng,mesh,DoF,FE);
% Allocate some terms for power iteration
% ------------------------------------------------------------------------------
keff = data.problem.KeffGuess;
oldkeff = 0.0; rhodenomphi = 1.0;
pwTolphi = data.solver.PIFluxTolerance;
pwTolkeff = data.solver.PIKeffTolerance;
% Perform power iteration
% ------------------------------------------------------------------------------
glob_fiss_rate = calculate_global_QoI(data.XS(xsid),mesh,DoF,FE,data.Fluxes.Phi,'Production');
for pi=1:data.solver.PIMaxIterations
    % Print PI Iteration Information
    if data.IO.PrintIterationInfo
        msg = sprintf('Begin Power Iteration #: %d',pi);
        disp(msg);
    end
    % Compute fission source
    data.Sources.FissionSource = SetFissionSource_Transport(data.XS(xsid),1:ng,1:ng,data.Fluxes.Phi,mesh,DoF,FE,keff);
    % Solve linear transport iteration
    data = solve_linear_iteration_transport(data,xsid,qid,mesh,DoF,FE);
    % Update fission source and keff
    old_glob_fiss_rate = glob_fiss_rate;
    glob_fiss_rate = calculate_global_QoI(data.XS(xsid),mesh,DoF,FE,data.Fluxes.Phi,'Production');
    oldoldkeff = oldkeff; oldkeff = keff;
    keff = oldkeff*glob_fiss_rate/old_glob_fiss_rate;
    % Compute spectral radius to guard for false convergence
    [rhonumerphi,pwErrphi] = compute_global_PI_error(data.Fluxes.Phi,KeffOldPhi);
    if data.solver.PISpectralRadiusCheck
        rhokeff = abs(keff-oldkeff)/abs(oldkeff-oldoldkeff);
        rhophi = rhonumerphi / rhodenomphi;
        rhodenomphi = rhonumerphi;
        if pi > 3
            pwTolkeff = (1.0-rhokeff)*data.solver.PIKeffTolerance;
            pwTolphi =  (1.0-rhophi) *data.solver.PIFluxTolerance;
        end
    end
    % Print PI iteration information
    pwErrkeff = abs(1.0 - oldkeff / keff);
    if data.IO.PrintIterationInfo
        msg = sprintf('  Power Iteration %d: (keff|Phi) Error - (%0.5e | %0.5e)',pi,pwErrkeff,pwErrphi);
        disp(msg);
    end
    % Check for convergence
    if pwErrkeff < pwTolkeff && pwErrphi < pwTolphi && pi > 4
        if data.IO.PrintIterationInfo
            msg = sprintf('End Power Iteration #: %d',pi);
            disp(msg);
        end
        break; 
    end
    % Update solution information
    KeffOldPhi = data.Fluxes.Phi;
end
% Process some iteration/output information
% ------------------------------------------------------------------------------
data.problem.KeffGuess = keff;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Iterative Procedure for k-eigenvalue diffusion problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = perform_diffusion_eignevalue_problem(data, mesh, DoF, FE)
% Get some additional information
% ------------------------------------------------------------------------------
xsid = data.Transport.XSID;
ng = data.Groups.NumberEnergyGroups;
% Set initial source terms
% ------------------------------------------------------------------------------
KeffOldPhi = data.Fluxes.Phi;

% Allocate some terms for power iteration
% ------------------------------------------------------------------------------
keff = data.problem.KeffGuess;
oldkeff = 0.0; oldoldkeff = 0.0;
pwErrphi = 0.0; pwErrkeff = 0.0;
pwTolphi = data.solver.PIFluxTolerance;
pwTolkeff = data.solver.PIKeffTolerance;
% Perform power iteration
% ------------------------------------------------------------------------------
glob_fiss_rate = calculate_global_QoI(data.XS(xsid),mesh,DoF,FE,data.Fluxes.Phi,'Production');
for pi=1:data.solver.PIMaxIterations
    % Compute fission source
    
    % Solve linear diffusion problem
    
    % Update fission source and keff
    
    % Check for convergence
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Miscellaneous Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_execution_function_handle(data)
if strcmpi(data.problem.TransportMethod,'transport')
    if strcmpi(data.problem.ProblemType,'sourcedriven')
        out = @perform_transport_source_driven_problem;
    elseif strcmpi(data.problem.ProblemType,'eigenvalue')
        out = @perform_transport_eignevalue_problem;
    end
elseif strcmpi(data.problem.TransportMethod,'diffusion')
    if strcmpi(data.problem.ProblemType,'sourcedriven')
        out = @perform_diffusion_source_driven_problem;
    elseif strcmpi(data.problem.ProblemType,'eigenvalue')
        out = @perform_diffusion_eignevalue_problem;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = createAMRMesh(data, geometry)
if data.AMR.RefineMesh && data.AMR.RefinementLevels > 0
    geometry = AMRGeometry(geometry);
    if isfield(data.problem,'AMRIrregularity')
        geometry.MaxIrregularity = data.AMR.AMRIrregularity;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = prepare_problem_execution(data, mesh, DoF)
if strcmpi(data.problem.TransportMethod,'transport')
    data = prepare_transport_execution(data, mesh, DoF);
elseif strcmpi(data.problem.TransportMethod,'diffusion')
    data = prepare_diffusion_execution(data, mesh, DoF);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = prepare_transport_execution(data, mesh, DoF)
% Determine Angle Sets
data = determine_angle_sets_Rev1(data, mesh);
% Determine Sweep Orderings
if data.Transport.PerformSweeps
    data = calculate_sweep_orderings(data, mesh);
end
% Allocate Solution Space
xsid = data.Transport.XSID; qid = data.Transport.QuadID;
data = solution_allocation(data, DoF.TotalDoFs, data.Groups.NumberEnergyGroups, data.Quadrature(qid).TotalFluxMoments);
data = set_incident_boundary_fluxes(data,xsid,qid,mesh,DoF);
if data.Transport.HasReflectingBoundary
    data = calculate_reflecting_angles_Rev1( data, mesh );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = prepare_diffusion_execution(data, mesh, DoF)
% Determine Angle Sets
data = determine_angle_sets_Rev1(data, mesh);
% Determine Sweep Orderings
if data.Transport.PerformSweeps
    data = calculate_sweep_orderings(data, mesh);
end
% Allocate Solution Space
data = solution_allocation(data, DoF.TotalDoFs, data.Groups.NumberEnergyGroups, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = solution_allocation(data, ndof, mf, nf)
data.Fluxes.Phi = cell(mf,nf);
% Set initial flux values
for i=1:mf
    for j=1:nf
        if isa(data.Fluxes.StartingSolution, 'char')
            if strcmp(data.Fluxes.StartingSolution, 'random')
                data.Fluxes.Phi{i,j} = rand(ndof,1);
            elseif strcmp(data.Fluxes.StartingSolution, 'zero')
                data.Fluxes.Phi{i,j} = zeros(ndof,1);
            elseif strcmp(data.Fluxes.StartingSolution, 'one')
                data.Fluxes.Phi{i,j} = ones(ndof,1);
            elseif strcmp(data.Fluxes.StartingSolution, 'function')
                data.Fluxes.Phi{i,j} = ndat.StartingSolutionFunction{i,j}(DoF.NodeLocations);
            else
                data.Fluxes.Phi{i,j} = rand(ndof,1);
            end
        elseif isa(data.Fluxes.StartingSolution, 'double')
            data.Fluxes.Phi{i,j} = ndat.StartingSolution*ones(ndof,1);
        else
            error('Cannot determine starting flux value type.');
        end
    end
end
% Set PhiOld
data.Fluxes.PhiOld = data.Fluxes.Phi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [abs_err,rel_err] = compute_global_PI_error(phi,phi0)
if size(phi,1) ~= size(phi0,1), error('Energy dimension of fluxes do not match.'); end
delt = 1e-12; ng = size(phi,1);
abs_err = 0; rel_err = 0; count = 0;
for g=1:ng
    rel_err = rel_err + sum(((phi{g,1} - phi0{g,1})./(phi{g,1}+delt)).^2);
    abs_err = abs_err + sum((phi{g,1} - phi0{g,1}).^2);
    count = count + length(phi{g,1});
end
abs_err = sqrt(abs_err/count);
rel_err = sqrt(rel_err/count);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mms, deg] = get_mms_information(data)
if data.MMS.PerformMMS
    mms = true;  deg = data.MMS.QuadOrder;
else
    mms = false; deg = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_problem_mesh(data, mesh, DoF, FE)
determine_refinement_cells(data, mesh, DoF, FE);
mesh.refine_mesh();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%