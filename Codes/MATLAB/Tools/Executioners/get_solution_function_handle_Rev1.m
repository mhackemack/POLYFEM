%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get Solution Function Handle
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
function out = get_solution_function_handle_Rev1(data)
ftype = data.problem.FEMType;
ttype = data.problem.TransportMethod;
if strcmpi(ttype, 'diffusion')
    if strcmpi(ftype, 'dfem')
        out = @perform_dfem_diffusion_Rev1;
    elseif strcmpi(ftype, 'cfem')
        out = @perform_cfem_diffusion_Rev1;
    end
elseif strcmpi(ttype, 'transport')
    if strcmpi(ftype, 'dfem')
        tt = data.Transport.TransportType;
        if strcmpi(tt, 'upwind')
            if data.Transport.PerformSweeps
                out = @exec_func_dfem_transport_sweep_Rev1;
            else
                out = @exec_func_dfem_transport_upwind_Rev1;
            end
        elseif strcmpi(tt, 'hybrid')
            out = @exec_func_dfem_transport_hybrid;
        end
    elseif strcmpi(ftype, 'cfem')
        out = @perform_cfem_transport;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%