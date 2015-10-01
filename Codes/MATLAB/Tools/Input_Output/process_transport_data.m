%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Transport Input Data
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
function data = process_transport_data(data)
% Process General Transport Data
% ------------------------------------------------------------------------------
if ~isfield(data,'Transport'), error('Transport field is required.'); end
if ~isfield(data.Transport,'XSID'), error('Need XS ID.'); end
if ~isfield(data.Transport,'TransportType'), error('Need a transport type.'); end
% Check input of specific transport types
if strcmpi(data.Transport.TransportType,'upwind')
    if ~isfield(data.Transport,'PerformSweeps')
        data.Transport.PerformSweeps = 0;
        data.Transport.VisualizeSweeping = 0;
    end
    if ~isfield(data.Transport,'VisualizeSweeping')
        data.Transport.VisualizeSweeping = 0;
    end
elseif strcmpi(data.Transport.TransportType,'hybrid')
    data.Transport.PerformSweeps = 0;
    data.Transport.VisualizeSweeping = 0;
    if ~isfield(data.Transport,'StabilizationMethod')
        data.Transport.StabilizationMethod = 'EGDG';
    end
    if ~isfield(data.Transport,'FluxStabilization')
        data.Transport.FluxStabilization = 2.0;
    end
    if ~isfield(data.Transport,'CurrentStabilization')
        data.Transport.CurrentStabilization = 1.0;
    end
end
% Process Angular Quadrature Data
% ------------------------------------------------------------------------------
if ~isfield(data,'Quadrature'), error('No angular quadrature field specified.'); end
if ~isfield(data.Transport,'QuadID'), error('Need main quadrature ID.'); end
data = process_angular_quadrature(data);
% Process Transport Acceleration Data
% ------------------------------------------------------------------------------
if isfield(data,'Acceleration')

end
% Process FEM Booleans
% ------------------------------------------------------------------------------
data.problem.FEMVolumeBools  = [1,1,0];
data.problem.FEMSurfaceBools = [1,0,0,0];
