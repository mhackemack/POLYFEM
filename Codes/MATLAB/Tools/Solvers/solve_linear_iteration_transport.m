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
function [data, sol] = solve_linear_iteration_transport(data,mesh,DoF,FE,sol)
% Get Iteration Information and Function Handles
% ------------------------------------------------------------------------------
num_gs = data.Neutronics.NumGroupSets;
inv_func = get_solution_function_handle(data);
ags_maxits = data.solver.AGSMaxIterations;
wgs_maxits = data.solver.WGSMaxIterations;
ags_rel_tol = data.solver.AGSRelativeTolerance;
wgs_rel_tol = data.solver.WGSRelativeTolerance;
ags_abs_tol = data.solver.AGSAbsoluteTolerance;
wgs_abs_tol = data.solver.WGSAbsoluteTolerance;
ags_gs_upscatter = data.Neutronics.GroupSetUpscattering;
wgs_accel_bools = data.Neutronics.Transport.WGSAccelerationBool;
ags_accel_bools = data.Neutronics.Transport.AGSAccelerationBool;
wgs_accel_resid = data.Neutronics.Transport.WGSAccelerationResidual;
ags_accel_resid = data.Neutronics.Transport.AGSAccelerationResidual;
% Build some data solution structures
% ------------------------------------------------------------------------------
DSA_Matrix = [];
AGS_OldPhi = [];
WGS_OldPhi = [];
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
        % Iterate within a group set
        for it=1:wgs_maxits(gs)
            % Perform within-group jacobi iteration
            [data.Neutronics, gsx] = inv_func(data.Neutronics,mesh,DoF,FE);
            % Perform Acceleration if necessary
            
            % Check convergence criteria
            
        end
        % Determine if current group set is now converged
        if ~ags_gs_upscatter(gs)
            if wgs_rel_converged && wgs_abs_converged, gs_converged(gs) = true; end
        end
    end
    % Perform Acceleration if necessary
    
    % Perform error tolerance checks
    
end
% 
% ------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxialary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,sol] = apply_transport_preconditioner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%