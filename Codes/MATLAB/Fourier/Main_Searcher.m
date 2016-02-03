%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Fourier Analysis Script - Search Method
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
data.Type = 'Search';
% outputs
data.Output.plotting_bool = true;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry_type = 'cart';
% x=[logspace(-3,0,55),logspace(0,2,141),logspace(2,3,45)];
% x=unique(x);
log_xmin = 0; log_xmax = 0; xnum = 1;
x = logspace(log_xmin, log_xmax, xnum);
dyz = [1];
% dyz = [1/100,1/64,1/16,1/4,1,4,16,64,100];
nx = 1;
ny = 1;
nz = 1;
% mat regions
mats(1).ID = 2;
mats(1).Region = [0,0;.5,0;.5,.5;0,.5];
mats(2).ID = 2;
mats(2).Region = [.5,.5;1,.5;1,1;.5,1];
% fem
data.Neutronics.TransportMethod = 'SI';
data.Neutronics.Transport.transportType = 'upwind';
data.problem.refineMesh = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.FEMLumping = false;
data.Neutronics.SpatialMethod = 'LAGRANGE';
data.Neutronics.FEMType = 'DFEM';
data.Neutronics.DSAType = 'MIP';
data.Neutronics.IP_Constant = 4;
% hybrid transport properties
data.Neutronics.Transport.StabilizationMethod = 'EGDG';
data.Neutronics.Transport.FluxStabilization = 2;
data.Neutronics.Transport.CurrentStabilization = 1;
% angular quadrature
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = [4];
data.Neutronics.Transport.PnOrder = 0;
% groups
data.Neutronics.numberEnergyGroups = 1;
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
data.NumberPhasePerDim = 5;
% end user input section
% ------------------------------------------------------------------------------
% Populate data and output structures
% -----------------------------------
print_FA_heading(data, x, dyz);
[data, inputs] = process_fourier_inputs( data, x, dyz, [nx, ny, nz], [] );
inputs = build_phase_transformation_matrix(data, inputs);
% Retrieve all spectrum data and postprocess
% ------------------------------------------
outputs = calculate_eigenspectrums(data, inputs);
process_fourier_outputs(data, inputs, outputs);

