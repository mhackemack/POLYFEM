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
function out = get_solution_function_handle(data)
ftype = data.Neutronics.FEMType;
ttype = data.Neutronics.transportMethod;
if strcmpi(ttype, 'diffusion')
    if strcmpi(ftype, 'dfem')
        out = @perform_dfem_diffusion;
    elseif strcmpi(ftype, 'cfem')
        out = @perform_cfem_diffusion;
    end
elseif strcmpi(ttype, 'transport')
    if strcmpi(ftype, 'dfem')
        tt = data.Neutronics.Transport.transportType;
        if strcmpi(tt, 'upwind')
            if data.Neutronics.Transport.performSweeps
                out = @exec_func_dfem_transport_sweep;
            else
                out = @exec_func_dfem_transport_upwind;
            end
        elseif strcmpi(tt, 'hybrid')
            out = @exec_func_dfem_transport_hybrid;
        end
    elseif strcmpi(ftype, 'cfem')
        out = @perform_cfem_transport;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%