%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve Linear Iteration
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = solve_linear_iteration_diffusion(data,mesh,DoF,FE)
% Get Iteration Information and Function Handles
% ------------------------------------------------------------------------------
num_gs = data.Groups.NumberGroupSets;
inv_func = get_solution_function_handle(data);
ags_maxits = data.solver.AGSMaxIterations;
wgs_maxits = data.solver.WGSMaxIterations;
ags_rel_tol = data.solver.AGSRelativeTolerance;
wgs_rel_tol = data.solver.WGSRelativeTolerance;
ags_abs_tol = data.solver.AGSAbsoluteTolerance;
wgs_abs_tol = data.solver.WGSAbsoluteTolerance;
% Build some data solution structures
% ------------------------------------------------------------------------------
diff_xsid   = 1; % always 1 for the diffusion cross sections
Diff_Matrix = [];
AGS_OldPhi  = [];
% Perform Scattering Kernel Convergence
% ------------------------------------------------------------------------------
gs_converged = false(num_gs, 1);
for m=1:ags_maxits
    % Loop through all group sets
    for gs=1:num_gs
        % Skip if group set is already converged
        if gs_converged(gs), continue; end
        % Set convergence booleans
        wgs_rel_converged = false;
        wgs_abs_converged = false;
        % Get some additional group set information
        grps = data.Groups.GroupSet{gs};
        % Iterate within a group set
        for it=1:wgs_maxits(gs)
            % Update Solution Vectors
            data = set_old_fluxes(data,grps);
            % Perform within-group jacobi iteration
            
            % Retreive convergence criteria
            [err_L2, norm_L2] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,data.Fluxes.PhiOld,grps,1,2);
            [err_inf, norm_inf] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,data.Fluxes.PhiOld,grps,1,inf);
            if err_L2 / norm_L2 < wgs_rel_tol(gs), wgs_rel_converged = true; else wgs_rel_converged = false; end
            if err_inf < wgs_abs_tol(gs), wgs_abs_converged = true; else wgs_abs_converged = false; end
            % Exit if convergence is met
            if wgs_rel_converged && wgs_abs_converged, break; end
        end
    end
    % Perform error tolerance checks
    [err_L2, norm_L2] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,AGS_OldPhi,grps,1,2);
    [err_inf, norm_inf] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,AGS_OldPhi,grps,1,inf);
    if err_L2 / norm_L2 < ags_rel_tol, ags_rel_converged = true; else ags_rel_converged = false; end
    if err_inf < ags_abs_tol, ags_abs_converged = true; else ags_abs_converged = false; end
    % Update Solution Vectors
    AGS_OldPhi = data.Fluxes.Phi;
    % Exit if convergence is met
    if ags_rel_converged && ags_abs_converged, break; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxialary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = set_old_fluxes(data,groups)
ng = length(groups);
% Loop through energy groups and moments to update fluxes
for g=1:ng
    data.Fluxes.PhiOld{groups(g)} = data.Fluxes.Phi{groups(g)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%