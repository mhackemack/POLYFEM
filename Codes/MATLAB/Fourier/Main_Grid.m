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
data.Output.plotting_bool = false;
data.Output.printing_bool = false;
data.Output.file_bool = false;
% geometry
data.problem.Dimension = 2;
data.geometry_type = 'cart';
% x=[logspace(-3,0,35),logspace(0,2,100),logspace(2,3,40)];
% x=unique(x);
% log_xmin = 0; log_xmax = 0; xnum = 1;
% x = logspace(log_xmin, log_xmax, xnum);
x = 2e2;
dyz = [1];
% dyz = [1/100,1/16,1/4,4,16,100];
nx = 2;
ny = 2;
nz = 1;
% mat regions
mats(1).ID = 2;
mats(1).Region = [0,0;.5,0;.5,.5;0,.5];
mats(2).ID = 2;
mats(2).Region = [.5,.5;1,.5;1,1;.5,1];
% fem
data.problem.refineMesh = false;
data.Neutronics.FEMDegree = 1;
data.Neutronics.FEMLumping = false;
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
% groups
data.Neutronics.numberEnergyGroups = 1;
% xs
c = 0.9999; sigt = [1.0;100.0];
% c = 0.9999; sigt = 1.0;
data.Neutronics.Transport.TotalXS = sigt;
data.Neutronics.Diffusion.DiffusionXS = (1/3)./sigt;
data.Neutronics.Transport.ScatteringXS = c*sigt;
data.Neutronics.Diffusion.AbsorbXS = (1-c)*sigt;
% bcs
data.Neutronics.Transport.BCFlags = glob.Periodic;
data.Neutronics.Transport.BCVals = 0.0;
% phase
data.NumberPhasePerDim = 51;
% end user input section
% ------------------------------------------------------------------------------
% Populate data and output structures
% -----------------------------------
print_FA_heading(data, x, dyz);
[data, inputs] = process_fourier_inputs( data, x, dyz, [nx, ny, nz], mats );
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
