%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Input Data
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
function [data, geometry] = process_user_input(data, geometry)
% Process Geometry
% ------------------------------------------------------------------------------
if isstruct(geometry) && isfield(geometry, 'geometry')
    geometry = geometry.geometry;
end
if geometry.Dimension ~= data.problem.Dimension
    error('Dimensionality does match between input and geometry.')
end
% This is possibly dangerous...
if data.problem.NumberMaterials ~= max(geometry.MatID)
    data.problem.NumberMaterials = max(geometry.MatID);
end
% Process XS Data
% ------------------------------------------------------------------------------
data.XS = process_xs_data(data.XS);
% Process General Neutronics Data
% ------------------------------------------------------------------------------
if strcmpi(data.problem.ProblemType, 'SourceDriven')
    
elseif strcmpi(data.problem.ProblemType, 'Eigenvalue')
    if ~strcmpi(data.Fluxes.StartingSolution, 'one')
        data.Fluxes.StartingSolution = 'one';
    end
    if ~isfield(data.problem,'KeffGuess')
        data.problem.KeffGuess = 1.0;
    end
end
% Process Transport Problem Data
% ------------------------------------------------------------------------------
if strcmpi(data.problem.TransportMethod, 'transport')
    if ~isfield(data,'Transport'), error('No Transport Field.'); end
    data = process_transport_data(data, geometry);
end
% Process Diffusion Problem Data
% ------------------------------------------------------------------------------
if strcmpi(data.problem.TransportMethod, 'diffusion')
    if ~isfield(data,'Diffusion'), error('No Diffusion Field.'); end
    data = process_diffusion_data(data, geometry);
end
% Process DoF/FE Data
% ------------------------------------------------------------------------------
if strcmpi(data.problem.SpatialMethod, 'ld')
    data.problem.DoFType = 1;
elseif strcmpi(data.problem.SpatialMethod, 'lagrange')
    data.problem.DoFType = 1;
else
    data.problem.DoFType = 2;
end
if ~isfield(data.problem, 'FEMLumping')
    data.problem.FEMLumping = false;
end
% Check AMR Data
% ------------------------------------------------------------------------------
data = process_amr_data(data);
% Check MMS Data
% ------------------------------------------------------------------------------
data = process_mms_data(data);
% ------------------------------------------------------------------------------