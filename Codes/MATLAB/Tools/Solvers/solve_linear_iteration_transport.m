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
function data = solve_linear_iteration_transport(data,trans_xsid,trans_qid,mesh,DoF,FE)
% Get Some GroupSet Information
% ------------------------------------------------------------------------------
num_gs = data.Groups.NumberGroupSets; 
all_groups = 1:data.Groups.NumberEnergyGroups;
tot_moments = data.Quadrature(trans_qid).TotalFluxMoments;
% Get Iteration Information and Function Handles
% ------------------------------------------------------------------------------
inv_func = get_solution_function_handle_Rev1(data);
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
ActiveAccelID   = 0;
ActiveAccelInfo = [];
DSA_Matrix      = [];
AGS_OldPhi      = data.Fluxes.Phi;
% Perform Scattering Kernel Convergence
% ------------------------------------------------------------------------------
gs_converged = false(num_gs, 1);
for m=1:ags_maxits
    % Print AGS Iteration Information
    if data.IO.PrintIterationInfo
        msg = sprintf('  Begin AGS Iteration #: %d',m);
        disp(msg);
    end
    % Loop through all group sets
    for gs=1:num_gs
        % Skip if group set is already converged
        if gs_converged(gs), continue; end
        % Print WGS Iteration Information
        if data.IO.PrintIterationInfo
            msg = sprintf('    Begin WGS Iteration for GS: %d',gs);
            disp(msg);
        end
        % Set WGS convergence booleans
        wgs_rel_converged = false;
        wgs_abs_converged = false;
        % Get some additional group set information
        wid  = wgs_accel_ids(gs);
        grps = data.Groups.GroupSets{gs};
        % Build Into-Group-Set Scattering
        in_grps = all_groups; in_grps(grps) = [];
        data.Sources.InScatteringSource = SetScatteringSource_Transport(data.XS(trans_xsid),data.Quadrature(trans_qid),grps,in_grps,data.Fluxes.Phi,mesh,DoF,FE);
        % Iterate within a group set
        for it=1:wgs_maxits(gs)
            % Update Solution Vectors
            data = set_old_fluxes(data,grps,1:data.Quadrature(trans_qid).TotalFluxMoments);
            gs_flux = cell(length(all_groups),tot_moments);
            gs_flux_old = cell(length(all_groups),tot_moments);
            % Build Within-Group-Set Scattering
            data.Sources.WithinScatteringSource = SetScatteringSource_Transport(data.XS(trans_xsid),data.Quadrature(trans_qid),grps,grps,data.Fluxes.Phi,mesh,DoF,FE);
            % Perform within-group jacobi iteration - no opp. reflecting boundaries
            % ---------------------------------------------------------------------
            if ~data.Transport.HasOpposingReflectingBoundary
                [data,tflux] = inv_func(data,trans_xsid,trans_qid,grps,mesh,DoF,FE);
                data = set_fluxes(data,grps,1:data.Quadrature(trans_qid).TotalFluxMoments,tflux);
            % Perform within-group jacobi iteration - w/ opp. reflecting boundaries
            % This converges the boundary angular fluxes for a given scattering source
            % ---------------------------------------------------------------------
            else
                it_abs_bool = false; it_rel_bool = false; it_count = 0;
                while ~it_abs_bool && ~it_rel_bool
                    [data,tflux] = inv_func(data,trans_xsid,trans_qid,grps,mesh,DoF,FE);
                    % 
                    for g=1:length(grps)
                        for mm=1:tot_moments
                            gs_flux{grps(g),mm} = tflux{g,mm};
                        end
                    end
                    %
                    [err_L2, norm_L2] = compute_flux_moment_differences(DoF,FE,gs_flux,gs_flux_old,grps,1,2);
                    [err_inf, norm_inf] = compute_flux_moment_differences(DoF,FE,gs_flux,gs_flux_old,grps,1,inf);
                    % 
                    for g=1:length(grps)
                        for mm=1:tot_moments
                            gs_flux_old{grps(g),mm} = gs_flux{grps(g),mm};
                        end
                    end
                    % Check convergence
                    if err_L2 / norm_L2 < 1e-15, it_rel_bool = true; else it_rel_bool = false; end
                    if err_inf < 1e-15, it_abs_bool = true; else it_abs_bool = false; end
                    it_count = it_count + 1;
                end
                data = set_fluxes(data,grps,1:data.Quadrature(trans_qid).TotalFluxMoments,tflux);
            end
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
                data = exec_func_RHS_DSA(data,ActiveAccelID,ActiveAccelInfo.XSID,mesh,DoF,FE,data.Fluxes.Phi,data.Fluxes.PhiOld);
%                 data = exec_func_RHS_DSA(data,ActiveAccelID,ActiveAccelInfo.XSID,mesh,DoF,FE);
                [data, DSA_Matrix] = perform_transport_acceleration(data,ActiveAccelID,mesh,DoF,FE,DSA_Matrix);
                data.Acceleration.Residual{ActiveAccelID} = [];
            end
            % Retreive convergence criteria
            [err_L2, norm_L2] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,data.Fluxes.PhiOld,grps,1,2);
            [err_inf, norm_inf] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,data.Fluxes.PhiOld,grps,1,inf);
            if err_L2 / norm_L2 < wgs_rel_tol(gs), wgs_rel_converged = true; else wgs_rel_converged = false; end
            if err_inf < wgs_abs_tol(gs), wgs_abs_converged = true; else wgs_abs_converged = false; end
            % Print WGS Iteration Information
            if data.IO.PrintIterationInfo
                msg = sprintf('      (%d) WGS Iteration %d: (L2|Linf) Error - (%0.5e | %0.5e)',gs,it,err_L2 / norm_L2,err_inf);
                disp(msg);
            end
            % Exit if convergence is met
            if wgs_rel_converged && wgs_abs_converged
                if data.IO.PrintIterationInfo
                    msg = sprintf('    End WGS Iteration for GS: %d',gs);
                    disp(msg);
                end
                break; 
            end
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
        data = exec_func_RHS_DSA(data,ActiveAccelID,ActiveAccelInfo.XSID,mesh,DoF,FE,data.Fluxes.Phi,AGS_OldPhi);
%         data = exec_func_RHS_DSA(data,ActiveAccelID,ActiveAccelInfo.XSID,mesh,DoF,FE);
        [data, DSA_Matrix] = perform_transport_acceleration(data,ActiveAccelID,mesh,DoF,FE,DSA_Matrix);
        data.Acceleration.Residual{ActiveAccelID} = [];
    end
    % Perform error tolerance checks
    [err_L2, norm_L2] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,AGS_OldPhi,1:data.Groups.NumberEnergyGroups,1,2);
    [err_inf, norm_inf] = compute_flux_moment_differences(DoF,FE,data.Fluxes.Phi,AGS_OldPhi,1:data.Groups.NumberEnergyGroups,1,inf);
    if err_L2 / norm_L2 < ags_rel_tol, ags_rel_converged = true; else ags_rel_converged = false; end
    if err_inf < ags_abs_tol, ags_abs_converged = true; else ags_abs_converged = false; end
    % Print AGS Iteration Information
    if data.IO.PrintIterationInfo
        msg = sprintf('  AGS Iteration %d: (L2|Linf) Error - (%0.5e | %0.5e)',m,err_L2 / norm_L2,err_inf);
        disp(msg);
    end
    % Update Solution Vectors
    AGS_OldPhi = data.Fluxes.Phi;
    % Exit if convergence is met
    if ags_rel_converged && ags_abs_converged
        if data.IO.PrintIterationInfo
            msg = sprintf('  End AGS Iteration #: %d',m);
            disp(msg);
        end
        break;
    end
end
% 
% ------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxialary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = set_fluxes(data, grps, mom, flux)
ng = length(grps);
nm = length(mom);
% Loop through energy groups and moments to update fluxes
for g=1:ng
    for m=1:nm
        data.Fluxes.Phi{grps(g),mom(m)} = flux{g,m};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = set_old_fluxes(data, grps, mom)
ng = length(grps);
nm = length(mom);
% Loop through energy groups and moments to update fluxes
for g=1:ng
    for m=1:nm
        data.Fluxes.PhiOld{grps(g),mom(m)} = data.Fluxes.Phi{grps(g),mom(m)};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%