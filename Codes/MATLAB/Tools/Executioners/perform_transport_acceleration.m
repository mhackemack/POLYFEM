%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Perform Transport Acceleration Step
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_transport_acceleration(data,accel_id,mesh,DoF,FE,src,A)
% Retrieve some data information
% ------------------------------
global glob
a_type = data.Acceleration.Info(accel_id).AccelerationType;
groups = data.Acceleration.Info(accel_id).Groups;
mom    = data.Acceleration.Info(accel_id).Moments;
xsid   = data.Acceleration.Info(accel_id).XSID;
[a_handle, is_dsa] = get_accel_function_handle(data.Acceleration.Info(accel_id));
% Perform DSA Preconditioning
if is_dsa
    % Get DSA system matrix if not set
    if isempty(A)
        A = a_handle(data,accel_id,xsid,mesh,DoF,FE);
    end
    % Compute error and apply correction based on DSA type
    dx = A\src;
    if a_type == glob.Accel_WGS_DSA || a_type == glob.Accel_AGS_TG
        evec = data.Acceleration.Info(accel_id).ErrorShape;
        % Loop through energy groups
        for g=1:length(groups)
            data.Fluxes.Phi{g,1} = data.Fluxes.Phi{g,1} + evec(g)*dx;
        end
    elseif a_type == glob.Accel_Fission_DSA
        
    end
% Perform TSA Preconditioning
else
    % nothing yet - this will be a relatively quick fix at a later date...
end
% Apply acceleration outputs
varargout{1} = data;
if is_dsa && nargout > 1
    varargout{2} = A;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxialiary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a_handle, is_dsa] = get_accel_function_handle(accel)
global glob
accel_type = accel.AccelerationType;
accel_disc = accel.DiscretizationType;
% Search through acceleration types
if      accel_type == glob.Accel_WGS_DSA || ...
        accel_type == glob.Accel_AGS_TG  || ...
        accel_type == glob.Accel_Fission_DSA
    is_dsa = true;
    % Switch between diffusion discretizations
    if accel_disc == glob.Accel_DSA_MIP
        a_handle = @exec_func_MIP_DSA;
    elseif accel_disc == glob.Accel_DSA_IP
        a_handle = @exec_func_IP_DSA;
    elseif accel_disc == glob.Accel_DSA_DCF
        a_handle = @exec_func_DCF_DSA;
    elseif accel_disc == glob.Accel_DSA_M4S
        a_handle = @exec_func_M4S_DSA;
    end
elseif  accel_type == glob.Accel_WGS_TSA || ...
        accel_type == glob.Accel_AGS_TTG
    is_dsa = false;
    if accel_disc == glob.Accel_TSA_DFEM
        a_handle = @exec_func_DFEM_TSA;
    elseif accel_disc == glob.Accel_TSA_CFEM
        a_handle = @exec_func_CFEM_TSA;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%