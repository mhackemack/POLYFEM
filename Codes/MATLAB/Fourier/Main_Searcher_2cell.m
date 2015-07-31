%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Fourier Analysis Script - Search Method (2 Cell)
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
glob = get_globals('Office');
glob.print_info = false;
% Define all user inputs
% ------------------------------------------------------------------------------
data.Type = 'Search';
% outputs
data.Output.plotting_bool = true;
data.Output.printing_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 1;
data.geometry_type = 'cart';
log_xmin = -3; log_xmax = 3; xnum = 150;
x = logspace(log_xmin, log_xmax, xnum);
dyz = [1];
% fem
data.problem.refineMesh = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.SpatialMethod = 'PWLD';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.TransportMethod = 'SI';
data.Neutronics.DSAType = 'MIP';
data.Neutronics.IP_Constant = 2;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [2];
data.Neutronics.Transport.fluxMoments = 0;
% xs
c = 0.9999; sigt1 = 1.0; sigt2 = 1.0;
data.Neutronics.Transport.TotalXS = [sigt1,sigt2];
data.Neutronics.Diffusion.DiffusionXS = [1/(3*sigt1), 1/(3*sigt2)];
data.Neutronics.Transport.ScatteringXS = c*[sigt1,sigt2];
data.Neutronics.Diffusion.AbsorbXS = (1-c)*[sigt1,sigt2];
% bcs
data.Neutronics.Transport.BCFlags = glob.Periodic;
data.Neutronics.Transport.BCVals = 0.0;
% phase
data.NumberPhasePerDim = 4;
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
process_fourier_outputs(data, inputs, outputs);

