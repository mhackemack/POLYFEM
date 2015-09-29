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
if strcmp(data.problem.ProblemType,'SourceDriven')
    pcall = @perform_source_driven_problem;
elseif strcmp(data.problem.ProblemType,'Eigenvalue')
    pcall = @perform_eignevalue_problem;
end
geometry = createAMRMesh(data, geometry);
[mms, deg] = get_mms_information(data);
% ------------------------------------------------------------------------------
% Perform Initial Solution Calculation (terminates if no mesh refinement)
% ------------------------------------------------------------------------------
% Get Degree of Freedom Handler - gets remade after a mesh refinement
DoF = DoFHandler(geometry, data.problem.FEMDegree, data.problem.FEMType, data.problem.DoFType);
% Get FEHandler based on MMS - this is necessary so that FEHandler
% generates Quadrature Sets for the Polygonal/Polyhedral cells
FE = FEHandler(geometry, DoF, data.problem.SpatialMethod, data.problem.FEMVolumeBools, data.problem.FEMSurfaceBools, mms, deg);
% Allocate Solution Space - gets reallocated after a mesh refinement
data = prepare_problem_execution(data, geometry);
[data, sol] = solution_allocation(data, geometry, DoF);
% Solve the neutronics problem (either source-driven or k-eigenvalue)
[data, sol] = pcall(data, geometry, DoF, FE, sol);
sol.CellVertexNumbers = geometry.CellVertexNumbers;
% Calculate MMS Error if necessary
if mms
    sol.MMS_error = calculate_MMS_error(data, geometry, DoF, FE, sol.flux);
end
% Plot solution for viewer to see
if data.IO.PlotSolution
    glob.plot_counter = glob.plot_counter + 1;
    figure(glob.plot_counter)
    plot_solution(geometry,DoF,FE,sol.flux)
end
% Save off output objects if necessary
if data.IO.SaveSolution
    f_name = ['outputs/',data.IO.Path,'/',data.IO.Name];
    o_str = '_0';
    save([f_name,'_data',o_str,'.mat'], 'data');
    save([f_name,'_geometry',o_str,'.mat'], 'geometry');
    save([f_name,'_DoF',o_str,'.mat'], 'DoF');
    save([f_name,'_FE',o_str,'.mat'], 'FE');
    save([f_name,'_sol',o_str,'.mat'], 'sol');
end
% Save off solution to VTK output for viewing in Visit/Paraview
if data.IO.SaveVTKSolution
    f_name = ['outputs/',data.IO.Path,'/',data.IO.Name];
    sol_flux = cell(data.Groups.NumberEnergyGroups,1);
    sol_name = cell(data.Groups.NumberEnergyGroups,1);
    for i=1:data.Groups.NumberEnergyGroups
        sol_flux{i} = sol.flux{i,1};
        sol_name{i} = ['flux_g',num2str(i)];
    end
    write_output_to_vtk_rev2(f_name, data, geometry, DoF, sol_flux, sol_name);
end
% ------------------------------------------------------------------------------
% Perform Mesh Refinement Calculations
% ------------------------------------------------------------------------------
if data.AMR.RefineMesh && data.AMR.RefinementLevels > 0
    % Retrieve some quick info first
    tsol = sol; clear sol;
    sol{1} = tsol;
    % Loop through refinement levels
    for iter=1:data.AMR.RefinementLevels
        o_str = ['_',num2str(iter)];
        r = iter + 1;
        % Save old objects for flux interpolation if necessary
        if data.AMR.ProjectSolution
            data0 = data;
            DoF0  = DoF;
        end
        % Form new mesh/DoFHandler/FEHandler
        refine_problem_mesh(data, geometry, DoF, FE, sol{r-1}.flux); % This works because of pass-by-reference
        DoF = DoFHandler(geometry, data.problem.FEMDegree, data.problem.FEMType, data.problem.DoFType);
        FE = FEHandler(geometry, DoF, data.problem.SpatialMethod, data.problem.FEMVolumeBools, data.problem.FEMSurfaceBools, mms, deg);
        data = prepare_problem_execution(data, geometry);
        [data, sol{r}] = solution_allocation(data, geometry, DoF);
        % Interpolate flux solutions to next refinement level
        if data.problem.projectSolution
            [data,sol{r}.flux] = interpolate_ref_flux(data0,data,geometry,DoF0,DoF,sol{r-1}.flux,sol{r}.flux);
        end
        % Compute new flux solution
        [data, tsol] = pcall(data, geometry, DoF, FE, sol{r});
        sol{r} = tsol;
        sol{r}.CellVertexNumbers = geometry.CellVertexNumbers;
        % Calculate MMS Error if necessary
        if mms
            sol{r}.MMS_error = calculate_MMS_error(data, geometry, DoF, FE, sol{r}.flux);
        end
        % Save off output objects if necessary
        if data.problem.SaveSolution
            save([f_name,'_data',o_str,'.mat'], 'data');
            save([f_name,'_geometry',o_str,'.mat'], 'geometry');
            save([f_name,'_DoF',o_str,'.mat'], 'DoF');
            save([f_name,'_FE',o_str,'.mat'], 'FE');
            % Solution Vector
            tsol = sol; sol = tsol{end};
            save([f_name,'_sol',o_str,'.mat'], 'sol');
            sol = tsol;
        end
        % Save off solution to VTK output for viewing in Visit/Paraview
        if data.problem.SaveVTKSolution
            f_name = ['outputs/',data.problem.Path,'/',data.problem.Name];
            sol_flux = cell(data.Neutronics.numberEnergyGroups,1);
            sol_name = cell(data.Neutronics.numberEnergyGroups,1);
            for i=1:data.Neutronics.numberEnergyGroups
                sol_flux{i} = sol{r}.flux{i,1};
                sol_name{i} = ['flux_g',num2str(i)];
            end
            write_output_to_vtk_rev2(f_name, data, geometry, DoF, sol_flux, sol_name);
        end
    end
end
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Iterative Procedure for source-driven neutronics problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, sol] = perform_source_driven_problem(data, geometry, DoF, FE, sol)
global glob
% Create DoF and other structures
% -------------------------------
ndat = data.Neutronics;
ndat.keff = 1.0;
solvdat = data.solver;
rT = solvdat.relativeTolerance;
aT = solvdat.absoluteTolerance;
fcall = get_solution_function_handle(data);
sol.flux0 = sol.flux;
mat = []; ferr0 = 1.0;

% Loop through iterations
% -----------------------
for l=1:solvdat.maxIterations
    if l > 1 && ~data.Neutronics.MultipleIterations, break; end
    if glob.print_info, disp(['-> Perform Source-Driven Iteration: ',num2str(l)]); end
    
    tictime = tic;
    [ndat, sol.flux, mat] = fcall(ndat, solvdat, geometry, DoF, FE, sol.flux, mat);
    [err_L2, norm_L2] = compute_flux_moment_differences(DoF, FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,2);
    [err_inf, norm_inf] = compute_flux_moment_differences(DoF, FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,inf);
    % Output iteration data
    if glob.print_info
        disp(['   -> Flux L2 Error: ',num2str(err_L2/norm_L2, '%0.9e')])
        disp(['   -> Flux L0 Error: ',num2str(err_inf, '%0.9e')])
        disp(' ')
    end
    % Check for Convergence
    toctime(l)   = toc(tictime);
    error_L2(l)  = err_L2;
    error_inf(l) = err_inf;
    n_L2(l)      = norm_L2;
    n_inf(l)     = norm_inf;
    if err_L2/norm_L2 < rT && err_inf < aT, break; end
    sol.flux0 = sol.flux;
end

% Apply outputs
% -------------
data.Neutronics = ndat;
sol.iter        = l;
sol.times       = toctime;
sol.error_L2    = error_L2;
sol.error_inf   = error_inf;
sol.norm_L2     = n_L2;
sol.norm_inf    = n_inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Iterative Procedure for k-eigenvalue neutronics problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = perform_eignevalue_problem(data, geometry, DoF, FE)
global glob
% Create DoF and other structures
% -------------------------------
ndat = data.Neutronics;
ndat.keff = 1.0;
solvdat = data.solver;
fcall = get_solution_function_handle(data);
rT = solvdat.relativeTolerance;
aT = solvdat.absoluteTolerance;

% Loop through iterations
% -----------------------
keff0 = ndat.keff;
mat = [];
sol.flux0 = sol.flux;
for l=1:solvdat.maxIterations
    if glob.print_info, disp(['-> Perform Eigenvalue Iteration: ',num2str(l)]); end
    
    tictime = tic;
    [ndat, sol.flux, mat] = fcall(ndat,solvdat,geometry,DoF,FE,sol.flux,mat);
    % Compute new keff and errors
    [keff,sol.flux] = estimate_new_keff(sol.flux,sol.flux0);
    [err_L2, norm_L2] = compute_flux_moment_differences(DoF,FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,2);
    [err_inf, norm_inf] = compute_flux_moment_differences(DoF, FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,inf);
    kerr(l) = abs(keff - keff0) / abs(keff);
    % Output iteration data
    if glob.print_info
        disp(['   -> keff:          ',num2str(keff)])
        disp(['   -> keff Error:    ',num2str(kerr(l), '%0.9e')])
        disp(['   -> Flux L2 Error: ',num2str(err_L2/norm_L2, '%0.9e')])
        disp(['   -> Flux L0 Error: ',num2str(err_inf, '%0.9e')])
        disp(' ')
    end
    % Check for Convergence
    toctime(l)   = toc(tictime);
    error_L2(l)  = err_L2;
    error_inf(l) = err_inf;
    n_L2(l)      = norm_L2;
    n_inf(l)     = norm_inf;
    if err_L2/norm_L2 < rT && err_inf< aT && kerr(l) < rT
        break
    end
    sol.flux0 = sol.flux;
    keff0 = keff;
    ndat.keff = keff;
    
end

% Apply outputs
% -------------
data.Neutronics = ndat;
sol.keff        = keff;
sol.iter        = l;
sol.times       = toctime;
sol.error_L2    = error_L2;
sol.error_inf   = error_inf;
sol.norm_L2     = n_L2;
sol.norm_inf    = n_inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Miscellaneous Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = createAMRMesh(data, geometry)
if data.AMR.RefineMesh && data.AMR.RefinementLevels > 0
    geometry = AMRGeometry(geometry);
    if isfield(data.problem,'AMRIrregularity')
        geometry.MaxIrregularity = data.AMR.AMRIrregularity;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = prepare_problem_execution(data, mesh)
% Prepares Angle Angregation and Sweep Orderings if necessary
if strcmp(data.problem.TransportMethod, 'Transport')
    % Reflecting Boundaries
    if data.Transport.HasReflectingBoundary
        data = determine_reflecting_boundaries( data, mesh );
    end
    % Determine Angle Sets
    data = determine_angle_sets(data, mesh);
    % Determine Sweep Orderings
    if data.Transport.PerformSweeps
        data = calculate_sweep_orderings(data, mesh);
    end
elseif strcmp(data.problem.TransportMethod, 'Diffusion')
    % nothing at this time...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = solution_allocation(data, mesh, DoF)
ndof = DoF.TotalDoFs;
% Generate DoF space for each energy group and flux moment
% --------------------------------------------------------
if strcmp(data.problem.TransportMethod, 'Diffusion')
    mf = data.Groups.NumberEnergyGroups;
    nf = 1;
else
    if data.Transport.HasReflectingBoundary
        data = allocate_reflecting_flux(data, mesh, DoF);
        data = calculate_reflecting_angles(data, mesh);
    end
    mf = data.Groups.NumberEnergyGroups;
    nf = data.Quadrature(1).TotalFluxMoments;
end
data.Fluxes.Phi = cell(mf,nf);
% Set initial flux values
for i=1:mf
    for j=1:nf
        if isa(ndat.StartingSolution, 'char')
            if strcmp(ndat.StartingSolution, 'random')
                data.Fluxes.Phi{i,j} = rand(ndof,1);
            elseif strcmp(ndat.StartingSolution, 'zero')
                data.Fluxes.Phi{i,j} = zeros(ndof,1);
            elseif strcmp(ndat.StartingSolution, 'one')
                data.Fluxes.Phi{i,j} = ones(ndof,1);
            elseif strcmp(ndat.StartingSolution, 'function')
                data.Fluxes.Phi{i,j} = ndat.StartingSolutionFunction{i,j}(DoF.NodeLocations);
            else
                data.Fluxes.Phi{i,j} = rand(ndof,1);
            end
        elseif isa(ndat.StartingSolution, 'double')
            data.Fluxes.Phi{i,j} = ndat.StartingSolution*ones(ndof,1);
        else
            error('Cannot determine starting flux value type.');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [keff,flux] = estimate_new_keff(flux,flux0)
ng = size(flux,1);
num = 0; denom = 0;
for g=1:ng
    num = num + norm(flux{g,1});
    denom = denom + norm(flux0{g,1});
end
keff = num/denom;
for g=1:ng
    flux{g,1} = flux{g,1} / keff;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = allocate_reflecting_flux(data, mesh, DoF)
global glob
% Set default ID information - this is true for head problems
xsid = 1;
qid = 1;
% Gather more data
mquad = data.Quadrature(qid);
na = mquad.NumberAngularDirections;
ng = data.Groups.NumberEnergyGroups;
% Clear then allocate reflecting fluxes (handles mesh refinement)
data.Fluxes.ReflectingFluxes = [];
data.Fluxes.ReflectingFluxes = cell(DoF.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fnorm = mesh.FaceNormal(ff,:);
    nflag = data.XS(xsid).BCFlags(mesh.FaceID(ff));
    if nflag == glob.Reflecting
        ndf = length(DoF.FaceCellNodes{ff,1});
        data.Fluxes.ReflectingFluxes{ff} = cell(na, 1);
        for m=1:mquad.NumberAngularDirections
            angDir = mquad.AngularDirections(m,:);
            if dot(fnorm, angDir) > 0
                data.Fluxes.ReflectingFluxes{ff}{m} = zeros(ndf, ng);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mms, deg] = get_mms_information(data)
if data.MMS.PerformMMS
    mms = true;
    deg = data.MMS.QuadOrder;
else
    mms = false;
    deg = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_problem_mesh(data, mesh, DoF, FE)
determine_refinement_cells(data, mesh, DoF, FE);
mesh.refine_mesh();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%