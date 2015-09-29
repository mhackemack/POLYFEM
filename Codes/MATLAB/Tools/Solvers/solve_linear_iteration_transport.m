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
function data = solve_linear_iteration_transport(data,mesh,DoF,FE)
% Get Iteration Information and Function Handles
% ------------------------------------------------------------------------------
num_gs = data.Groups.NumberGroupSets;
[inv_func, rhs_func] = get_solution_function_handle(data);
ags_maxits = data.solver.AGSMaxIterations;
wgs_maxits = data.solver.WGSMaxIterations;
ags_rel_tol = data.solver.AGSRelativeTolerance;
wgs_rel_tol = data.solver.WGSRelativeTolerance;
ags_abs_tol = data.solver.AGSAbsoluteTolerance;
wgs_abs_tol = data.solver.WGSAbsoluteTolerance;
ags_gs_upscatter = data.Groups.GroupSetUpscattering;
wgs_accel_bools = data.Acceleration.WGSAccelerationBool;
ags_accel_bool  = data.Acceleration.AGSAccelerationBool;
wgs_accel_resid = data.Acceleration.WGSAccelerationResidual;
ags_accel_resid = data.Acceleration.AGSAccelerationResidual;
wgs_accel_ids   = data.Acceleration.WGSAccelerationID;
ags_accel_id    = data.Acceleration.AGSAccelerationID;
% Build some data solution structures
% ------------------------------------------------------------------------------
trans_xsid = 1; % always 1 for the Transport cross sections
trans_qid = 1; % always 1 for the Transport Quadrature Set
ActiveAccelID   = 0;
ActiveAccelInfo = [];
DSA_Matrix      = [];
AGS_OldPhi      = [];
% Perform Scattering Kernel Convergence
% ------------------------------------------------------------------------------
gs_converged = false(num_gs, 1);
for m=1:ags_maxits
    % Loop through all group sets
    for gs=1:num_gs
        % Skip if group set is already converged
        if gs_converged(gs), continue; end
        % Set WGS convergence booleans
        wgs_rel_converged = false;
        wgs_abs_converged = false;
        % Get some additional group set information
        wid  = wgs_accel_ids(gs);
        grps = data.Groups.GroupSet{gs};
        % Iterate within a group set
        for it=1:wgs_maxits(gs)
            % Perform within-group jacobi iteration
            src = rhs_func(data,trans_xsid,trans_qid,grps,mesh,DoF,FE);
            data = inv_func(data,trans_xsid,trans_qid,grps,mesh,DoF,FE,src);
            % Gather Acceleration Residual if necessary
            if wgs_accel_resid(gs)
                
            end
            % Perform Acceleration if necessary
            if wgs_accel_bools(gs)
                % Update acceleration info if necessary
                if ActiveAccelID ~= wid
                    DSA_Matrix = [];
                    ActiveAccelID = wid;
                    ActiveAccelInfo = data.Acceleration.Info(wid);
                end
                % Perform Acceleration Here
                data = exec_func_RHS_DSA(data,ActiveAccelID,ActiveAccelInfo.XSID,mesh,DoF,FE);
                [data, DSA_Matrix] = perform_transport_acceleration(data,ActiveAccelID,mesh,DoF,FE,DSA_Matrix);
                data.Acceleration.Residual{ActiveAccelID} = [];
            end
            % Retreive convergence criteria
            [err_L2, norm_L2] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,data.Fluxes.PhiOld,grps,1,2);
            [err_inf, norm_inf] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,data.Fluxes.PhiOld,grps,1,inf);
            if err_L2 / norm_L2 < wgs_rel_tol(gs), wgs_rel_converged = true; else wgs_rel_converged = false; end
            if err_inf < wgs_abs_tol(gs), wgs_abs_converged = true; else wgs_abs_converged = false; end
            % Update Solution Vectors
            data = set_old_fluxes(data,grps,1:data.Quadrature(trans_qid).TotalFluxMoments);
            % Exit if convergence is met
            if wgs_rel_converged && wgs_abs_converged, break; end
        end
        % Determine if current group set is now converged
        if ~ags_gs_upscatter(gs)
            if wgs_rel_converged && wgs_abs_converged, gs_converged(gs) = true; end
        end
    end
    % Gather Acceleration Residual if necessary
    if ags_accel_resid
        
    end
    % Perform Acceleration if necessary
    if ags_accel_bool
        % Update acceleration info if necessary
        if ActiveAccelID ~= ags_accel_id
            DSA_Matrix = [];
            ActiveAccelID = ags_accel_id;
            ActiveAccelInfo = data.Acceleration.Info(ags_accel_id);
        end
        % Perform Acceleration Here
        data = exec_func_RHS_DSA(data,ActiveAccelID,ActiveAccelInfo.XSID,mesh,DoF,FE);
        [data, DSA_Matrix] = perform_transport_acceleration(data,ActiveAccelID,mesh,DoF,FE,DSA_Matrix);
        data.Acceleration.Residual{ActiveAccelID} = [];
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
% 
% ------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxialary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = set_old_fluxes(data,groups,mom)
ng = length(groups);
nm = length(mom);
% Loop through energy groups and moments to update fluxes
for g=1:ng
    for m=1:nm
        data.Fluxes.PhiOld{g,m} = data.Fluxes.Phi{g,m};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%