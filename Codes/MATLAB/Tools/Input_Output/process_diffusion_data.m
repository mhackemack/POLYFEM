%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Diffusion Input Data
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
function data = process_diffusion_data(data)
% Process General Diffusion Data
% ------------------------------------------------------------------------------
if ~isfield(data,'Diffusion'), error('Diffusion field is required.'); end

% Process FEM Booleans
% ------------------------------------------------------------------------------
data.problem.FEMVolumeBools  = logical([1,0,1]);
data.problem.FEMSurfaceBools = logical([1,0,0,0]);
if strcmpi(data.problem.FEMType,'dfem'), data.problem.FEMSurfaceBools(2) = true; end