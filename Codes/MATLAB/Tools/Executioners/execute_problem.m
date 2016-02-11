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
function [data, sol, geometry, DoF, FE] = execute_problem(data, geometry)
global glob
% Generate problem execution schedule
% ------------------------------------------------------------------------------
pcall = get_problem_function_handles(data);
geometry = createAMRMesh(data, geometry);
[mms, deg] = get_mms_information(data);
% ------------------------------------------------------------------------------
% Perform Initial Solution Calculation (terminates if no mesh refinement)
% ------------------------------------------------------------------------------
% Get Degree of Freedom Handler - gets remade after a mesh refinement
DoF = DoFHandler(geometry, data.Neutronics.FEMDegree, data.Neutronics.FEMType, data.Neutronics.DoFType);
% Get FEHandler based on MMS - this is necessary so that FEHandler
% generates Quadrature Sets for the Polygonal/Polyhedral cells
FE = FEHandler(geometry, DoF, data.Neutronics.SpatialMethod, data.Neutronics.FEMLumping, data.Neutronics.FEMVolumeBools, data.Neutronics.FEMSurfaceBools, mms, deg);
% Allocate Solution Space - gets reallocated after a mesh refinement
data = prepare_problem_execution(data, geometry);
[data.Neutronics, sol] = solution_allocation(data.Neutronics, geometry, DoF);
% Solve the neutronics problem (either source-driven or k-eigenvalue)
[data, sol] = pcall(data, geometry, DoF, FE, sol);
sol.CellVertexNumbers = geometry.CellVertexNumbers;
% Compute average material fluxes and some qoi's
sol.AverageMaterialFlux = calculate_average_material_QoI(data,geometry,DoF,FE,sol.flux);
sol.TotalMaterialFlux = calculate_total_QoI(data,geometry,DoF,FE,sol.flux,'Flux');
sol.TotalMaterialInteraction = calculate_total_QoI(data,geometry,DoF,FE,sol.flux,'Total');
sol.TotalMaterialAbsorption = calculate_total_QoI(data,geometry,DoF,FE,sol.flux,'Absorption');
% Calculate MMS Error if necessary
if mms, sol.MMS_error = calculate_MMS_error(data, geometry, DoF, FE, sol.flux); end
% Plot solution for viewer to see
if data.problem.plotSolution
    glob.plot_counter = glob.plot_counter + 1;
    figure(glob.plot_counter)
    plot_solution(geometry,DoF,FE,sol.flux)
end
% Save off output objects if necessary
if data.problem.saveSolution
    f_name = ['outputs/',data.problem.Path,'/',data.problem.Name];
    o_str = '_0';
    save([f_name,'_data',o_str,'.mat'], 'data');
    save([f_name,'_geometry',o_str,'.mat'], 'geometry');
    save([f_name,'_DoF',o_str,'.mat'], 'DoF');
    save([f_name,'_FE',o_str,'.mat'], 'FE');
    save([f_name,'_sol',o_str,'.mat'], 'sol');
end
% Save off solution to VTK output for viewing in Visit/Paraview
if data.problem.saveVTKSolution
    f_name = ['outputs/',data.problem.Path,'/',data.problem.Name];
    sol_flux = cell(data.Neutronics.numberEnergyGroups,1);
    sol_name = cell(data.Neutronics.numberEnergyGroups,1);
    for i=1:data.Neutronics.numberEnergyGroups
        sol_flux{i} = sol.flux{i,1};
        sol_name{i} = ['flux_g',num2str(i)];
    end
    write_output_to_vtk_rev2(f_name, data, geometry, DoF, sol_flux, sol_name);
end
% ------------------------------------------------------------------------------
% Perform Mesh Refinement Calculations
% ------------------------------------------------------------------------------
if data.problem.refineMesh && data.problem.refinementLevels > 0
    % Retrieve some quick info first
    tsol = sol; clear sol;
    sol{1} = tsol;
    % Loop through refinement levels
    for iter=1:data.problem.refinementLevels
        o_str = ['_',num2str(iter)];
        r = iter + 1;
        % Save old objects for flux interpolation if necessary
        if data.problem.projectSolution
            DoF0  = DoF;
            FE0 = FE;
        end
        % Form new mesh/DoFHandler/FEHandler
        refine_problem_mesh(data, geometry, DoF, FE, sol{r-1}.flux); % This works because of pass-by-reference
        DoF = DoFHandler(geometry, data.Neutronics.FEMDegree, data.Neutronics.FEMType, data.Neutronics.DoFType);
        FE = FEHandler(geometry, DoF, data.Neutronics.SpatialMethod, data.Neutronics.FEMLumping, data.Neutronics.FEMVolumeBools, data.Neutronics.FEMSurfaceBools, mms, deg);
        data = prepare_problem_execution(data, geometry);
        [data.Neutronics, sol{r}] = solution_allocation(data.Neutronics, geometry, DoF);
        % Interpolate flux solutions to next refinement level
        if data.problem.projectSolution
            sol{r}.flux = interpolate_ref_flux_Rev1(geometry,DoF0,DoF,FE0,sol{r-1}.flux);
        end
        % Compute new flux solution
        [data, tsol] = pcall(data, geometry, DoF, FE, sol{r});
        sol{r} = tsol;
        sol{r}.CellVertexNumbers = geometry.CellVertexNumbers;
        % Compute average material fluxes and some qoi's
        sol{r}.AverageMaterialFlux = calculate_average_material_QoI(data,geometry,DoF,FE,sol{r}.flux);
        sol{r}.TotalMaterialFlux = calculate_total_QoI(data,geometry,DoF,FE,sol{r}.flux,'Flux');
        sol{r}.TotalMaterialInteraction = calculate_total_QoI(data,geometry,DoF,FE,sol{r}.flux,'Total');
        sol{r}.TotalMaterialAbsorption = calculate_total_QoI(data,geometry,DoF,FE,sol{r}.flux,'Absorption');
        % Calculate MMS Error if necessary
        if mms, sol{r}.MMS_error = calculate_MMS_error(data, geometry, DoF, FE, sol{r}.flux); end
        % Save off output objects if necessary
        if data.problem.saveSolution
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
        if data.problem.saveVTKSolution
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
% ------------------------------------------------------------------------------
ndat = data.Neutronics;
ndat.keff = 1.0;
solvdat = data.solver;
rT = solvdat.relativeTolerance;
aT = solvdat.absoluteTolerance;
fcall = get_solution_function_handle(data);
sol.flux0 = sol.flux;
mat = []; ferr0 = 1.0;
% Loop through iterations
% ------------------------------------------------------------------------------
for l=1:solvdat.maxIterations
    if l > 1 && ~data.Neutronics.MultipleIterations, break; end
    if glob.print_info, disp(['-> Perform Source-Driven Iteration: ',num2str(l)]); end
    
    tictime = tic;
    [ndat, sol.flux, mat, extra_its, extra_time] = fcall(ndat, solvdat, geometry, DoF, FE, sol.flux, mat);
    [err_L2, norm_L2] = compute_flux_moment_differences(DoF, FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,2);
    [err_inf, norm_inf] = compute_flux_moment_differences(DoF, FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,inf);
%     [err_pw, norm_pw] = compute_flux_moment_differences(DoF, FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,'pdt_pw');
    % Output iteration data
    if glob.print_info
        disp(['   -> Flux L2 Error: ',num2str(err_L2/norm_L2, '%0.9e')])
        disp(['   -> Flux L0 Error: ',num2str(err_inf, '%0.9e')])
        disp(' ')
    end
    % Check for Convergence
    toctime(l)   = toc(tictime);
    DSA_its(l)   = extra_its;
    DSA_times(l) = extra_time;
    error_L2(l)  = err_L2;
    error_inf(l) = err_inf;
    n_L2(l)      = norm_L2;
    n_inf(l)     = norm_inf;
    if err_L2/norm_L2 < rT && err_inf< aT, break; end
    sol.flux0 = sol.flux;
end
% Apply outputs
% ------------------------------------------------------------------------------
data.Neutronics = ndat;
sol.iter        = l;
sol.DSA_iters   = DSA_its;
sol.times       = toctime;
sol.error_L2    = error_L2;
sol.error_inf   = error_inf;
sol.norm_L2     = n_L2;
sol.norm_inf    = n_inf;
% ------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Iterative Procedure for k-eigenvalue neutronics problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, sol] = perform_eignevalue_problem(data, geometry, DoF, FE, sol)
global glob
% Create DoF and other structures
% ------------------------------------------------------------------------------
ndat = data.Neutronics;
ndat.keff = 1.0;
solvdat = data.solver;
fcall = get_solution_function_handle(data);
rT = solvdat.relativeTolerance;
aT = solvdat.absoluteTolerance;
% Loop through iterations
% ------------------------------------------------------------------------------
keff0 = ndat.keff;
mat = [];
sol.flux0 = sol.flux;
for l=1:solvdat.maxIterations
    if glob.print_info, disp(['-> Perform Eigenvalue Iteration: ',num2str(l)]); end
    
    tictime = tic;
    [ndat, sol.flux, mat, extra_its, extra_time] = fcall(ndat,solvdat,geometry,DoF,FE,sol.flux,mat);
    % Compute new keff and errors
    [keff,sol.flux] = estimate_new_keff(sol.flux,sol.flux0);
    [err_L2, norm_L2] = compute_flux_moment_differences(DoF, FE,sol.flux,sol.flux0,1:sol.numberEnergyGroups,1,2);
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
    DSA_its(l)   = extra_its;
    DSA_times(l) = extra_time;
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
% ------------------------------------------------------------------------------
data.Neutronics = ndat;
sol.keff        = keff;
sol.iter        = l;
sol.DSA_iters   = DSA_its;
sol.times       = toctime;
sol.error_L2    = error_L2;
sol.error_inf   = error_inf;
sol.norm_L2     = n_L2;
sol.norm_inf    = n_inf;
% ------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Miscellaneous Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pcall = get_problem_function_handles(data)
if strcmp(data.problem.problemType,'SourceDriven')
    pcall = @perform_source_driven_problem;
elseif strcmp(data.problem.problemType,'Eigenvalue')
    pcall = @perform_eignevalue_problem;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = createAMRMesh(data, geometry)
if data.problem.refineMesh && data.problem.refinementLevels > 0
    geometry = AMRGeometry(geometry);
    if isfield(data.problem,'AMRIrregularity')
        geometry.MaxIrregularity = data.problem.AMRIrregularity;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = prepare_problem_execution(data, mesh)
% Prepares Angle Angregation and Sweep Orderings if necessary
if strcmp(data.Neutronics.transportMethod, 'Transport')
    if data.Neutronics.Transport.HasReflectingBoundary
        data = determine_reflecting_boundaries( data, mesh );
    end
    data = determine_angle_sets(data, mesh);
    if data.Neutronics.Transport.performSweeps
        data = calculate_sweep_orderings(data, mesh);
    end
elseif strcmp(data.Neutronics.transportMethod, 'Diffusion')
    % nothing at this time...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ndat, sol] = solution_allocation(ndat, mesh, DoF)
ndof = DoF.TotalDoFs;
% Generate DoF space for each energy group and flux moment
% --------------------------------------------------------
if strcmp(ndat.transportMethod, 'Diffusion')
    if ndat.Diffusion.HasPeriodicBoundary
        ndat = allocate_periodic_flux(ndat, mesh, DoF);
    end
    mf = ndat.numberEnergyGroups;
    nf = 1;
    ndat.TotalFluxMoments = 1;
else
    if ndat.Transport.HasReflectingBoundary
        ndat = allocate_reflecting_flux(ndat, mesh, DoF);
        ndat = calculate_reflecting_angles(ndat, mesh);
    end
    if ndat.Transport.HasPeriodicBoundary
        ndat = allocate_periodic_flux(ndat, mesh, DoF);
    end
    ndat = allocate_boundary_currents(ndat, mesh, DoF);
    ndat = allocate_beam_boundary_condition( ndat, mesh );
    mf = ndat.numberEnergyGroups;
    nf = ndat.TotalFluxMoments;
end
sol.flux = cell(mf,nf);
% Set initial flux values
for i=1:mf
    for j=1:nf
        if isa(ndat.StartingSolution, 'char')
            if strcmp(ndat.StartingSolution, 'random')
                sol.flux{i,j} = rand(ndof,1);
            elseif strcmp(ndat.StartingSolution, 'zero')
                sol.flux{i,j} = zeros(ndof,1);
            elseif strcmp(ndat.StartingSolution, 'one')
                sol.flux{i,j} = ones(ndof,1);
            elseif strcmp(ndat.StartingSolution, 'function')
                sol.flux{i,j} = ndat.StartingSolutionFunction{i,j}(DoF.NodeLocations);
            else
                sol.flux{i,j} = rand(ndof,1);
            end
        elseif isa(ndat.StartingSolution, 'double')
            sol.flux{i,j} = ndat.StartingSolution*ones(ndof,1);
        else
            error('Cannot determine starting flux value type.');
        end
    end
end
sol.numberEnergyGroups = mf;
sol.TotalFluxMoments = nf;
sol.SpatialDoFs = ndof;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, denom] = compute_flux_moment_differences(DoF, FE, flux, flux0, ngs, mom, n_type)
% global glob
err = 0; denom = 0;
if isinf(n_type)
    for g=1:length(ngs)
        gg = ngs(g);
        for m=1:length(mom)
            terr = max(abs(flux{gg,m} - flux0{gg,m}));
            tden = max(abs(flux{gg,m}));
            if terr > err, err = terr; end
            if tden > denom, denom = tden; end
        end
    end
elseif isa(n_type, 'double') && n_type > 0
    % Loop through Cells and compute L2 Norm
    for c=1:DoF.TotalCells
        M = FE.CellMassMatrix{c};
        cdofs = DoF.ConnectivityArray{c};
        zn = ones(length(cdofs), 1);
        for g=1:length(ngs)
            gg = ngs(g);
            for m=1:length(mom)
                err = err + (M*(flux{gg,m}(cdofs) - flux0{gg,m}(cdofs)).^n_type)'*zn;
                denom = denom + (M*flux{gg,m}(cdofs).^n_type)'*zn;
            end
        end
    end
    err = sqrt(err);
    denom = sqrt(denom);
elseif isa(n_type, 'char')
    if strcmpi(n_type, 'pdt_pw')
        for g=1:length(ngs)
            gg = ngs(g);
            sflux_max = max(abs([flux{gg,1}';flux0{gg,1}']))';
            for m=1:length(mom)
                terr = max(abs(flux{gg,m} - flux0{gg,m})./sflux_max);
                if terr > err, err = terr; end
%                 if tden > denom, denom = tden; end
            end
        end
    end
else
    error('Cannot determine norm type.');
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
function ndat = allocate_reflecting_flux(ndat, mesh, DoF)
global glob
na = ndat.Transport.NumberAngularDirections;
ng = ndat.numberEnergyGroups;
% Clear then allocate reflecting fluxes (handles mesh refinement)
ndat.Transport.ReflectingFluxes = [];
ndat.Transport.ReflectingFluxes = cell(DoF.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fnorm = mesh.FaceNormal(ff,:);
    nflag = ndat.Transport.BCFlags(mesh.FaceID(ff));
    if nflag == glob.Reflecting
        ndf = length(DoF.FaceCellNodes{ff,1});
        ndat.Transport.ReflectingFluxes{ff} = cell(na, 1);
        for m=1:ndat.Transport.NumberAngularDirections
            angDir = ndat.Transport.AngularDirections(m,:);
            if dot(fnorm, angDir) > 0
                ndat.Transport.ReflectingFluxes{ff}{m} = zeros(ndf, ng);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = allocate_periodic_flux(ndat, mesh, DoF)
global glob
ng = ndat.numberEnergyGroups;
if strcmp(ndat.transportMethod, 'Diffusion')
    ndat.Diffusion.PeriodicFluxes = [];
    ndat.Transport.PeriodicFluxes = cell(DoF.TotalFaces, 1);
    for f=1:mesh.TotalBoundaryFaces
        ff = mesh.BoundaryFaces(f);
        nfd = length(DoF.FaceCellNodes{ff,1});
        ndat.Transport.PeriodicFluxes{ff} = zeros(nfd, ng);
    end
elseif strcmp(ndat.transportMethod, 'Transport')
    na = ndat.Transport.NumberAngularDirections;
    ndat.Transport.PeriodicFluxes = [];
    ndat.Transport.PeriodicFluxes = cell(DoF.TotalFaces, 1);
    for f=1:mesh.TotalBoundaryFaces
        ff = mesh.BoundaryFaces(f);
        fnorm = mesh.FaceNormal(ff,:);
        nflag = ndat.Transport.BCFlags(mesh.FaceID(ff));
        if nflag == glob.Periodic
            ndf = length(DoF.FaceCellNodes{ff,1});
            ndat.Transport.PeriodicFluxes{ff} = cell(na, 1);
            for m=1:ndat.Transport.NumberAngularDirections
                angDir = ndat.Transport.AngularDirections(m,:);
                if dot(fnorm, angDir) < 0
                    ndat.Transport.PeriodicFluxes{ff}{m} = zeros(ndf, ng);
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = allocate_boundary_currents(ndat, mesh, DoF)
ng = ndat.numberEnergyGroups;
ndat.Transport.OutgoingCurrents = cell(DoF.TotalFaces, 1);
ndat.Transport.IncomingCurrents = cell(DoF.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    ndf = length(DoF.FaceCellNodes{ff,1});
    ndat.Transport.OutgoingCurrents{ff} = zeros(ndf, ng);
    ndat.Transport.IncomingCurrents{ff} = zeros(ndf, ng);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = allocate_beam_boundary_condition( ndat, mesh )
global glob
na = ndat.Transport.NumberAngularDirections;
ng = ndat.numberEnergyGroups;
ndat.Transport.BeamFluxes = cell(mesh.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fid = mesh.FaceID(ff);
    nflag = ndat.Transport.BCFlags(fid);
    if nflag == glob.IncidentBeam
        fnorm = mesh.FaceNormal(ff,:);
        ndat.Transport.BeamFluxes{ff} = zeros(na, ng);
        % Determine Maximum Angle
        min_dot = 0.0; min_ind = 1;
        for m=1:ndat.Transport.NumberAngularDirections
            angDir = ndat.Transport.AngularDirections(m,:);
            fdot = dot(fnorm, angDir);
            if fdot < 0 && fdot < min_dot && abs(fdot - min_dot) > 1e-12
                min_dot = fdot; min_ind = m;
            end
        end
        for g=1:ng
            ndat.Transport.BeamFluxes{ff}(min_ind,g) = ndat.Transport.BCVals{fid,g};
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mms, deg] = get_mms_information(data)
ttype = data.Neutronics.transportMethod;
if strcmp(ttype, 'Diffusion')
    if data.Neutronics.Diffusion.MMS
        mms = true;
        deg = data.Neutronics.Diffusion.QuadOrder;
    else
        mms = false;
        deg = 0;
    end
elseif strcmp(ttype, 'Transport')
    if data.Neutronics.Transport.MMS
        mms = true;
        deg = data.Neutronics.Transport.QuadOrder;
    else
        mms = false;
        deg = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_problem_mesh(data, mesh, DoF, FE, flux)
determine_refinement_cells(data, mesh, DoF, FE, flux);
mesh.refine_mesh();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%