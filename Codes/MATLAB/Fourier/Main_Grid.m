%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Fourier Analysis Script
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
% Clear Project Space
% -------------------
if exist('pbool', 'var')
    clearvars -except pbool
else
    clear; pbool = false;
end
clc; close all; %format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Define Path
% -----------
global glob
glob = get_globals('Home');
glob.print_info = false;
% Define all user inputs
% ------------------------------------------------------------------------------
data.Type = 'Grid';
% outputs
data.Output.plotting_bool = true;
data.Output.printing_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry_type = 'tri';
% x=[logspace(-3,0,35),logspace(0,2,100),logspace(2,3,40)];
% x=unique(x);
log_xmin = -1; log_xmax = -1; xnum = 1;
x = logspace(log_xmin, log_xmax, xnum);
dyz = [1];
% dyz = [1/100,1/16,1/4,4,16,100];
% fem
data.problem.refineMesh = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.TransportMethod = 'SI';
data.Neutronics.Transport.transportType = 'upwind';
data.Neutronics.DSAType = 'MIP';
data.Neutronics.IP_Constant = 4;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [4];
data.Neutronics.Transport.PnOrder = 0;
% xs
c = 0.9999; sigt = 1.0;
data.Neutronics.Transport.TotalXS = sigt;
data.Neutronics.Diffusion.DiffusionXS = 1/(3*sigt);
data.Neutronics.Transport.ScatteringXS = c*sigt;
data.Neutronics.Diffusion.AbsorbXS = (1-c)*sigt;
% bcs
data.Neutronics.Transport.BCFlags = glob.Periodic;
data.Neutronics.Transport.BCVals = 0.0;
% phase
data.NumberPhasePerDim = 151;
% end user input section
% ------------------------------------------------------------------------------
% Populate data and output structures
% -----------------------------------
print_FA_heading(data, x, dyz);
[data, inputs] = process_fourier_inputs( data, x, dyz );
inputs = build_phase_transformation_matrix(data, inputs);
% Retrieve all spectrum data and postprocess
% ------------------------------------------
outputs = calculate_eigenspectrums(data, inputs);
% Loop through quadrature sets
for q=1:length(data.Neutronics.Transport.SnLevels)
    % Loop through meshes
    for i=1:inputs.TotalMeshes
        
    end
end
