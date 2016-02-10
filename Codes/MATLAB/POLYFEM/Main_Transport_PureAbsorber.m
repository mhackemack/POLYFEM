%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D Pure Absorber Problem Run Script
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
% Prepare Project Space
% ------------------------------------------------------------------------------
clc; close all; format long e
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Office');
inp = 'Transport_PureAbsorber';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
geom_in.Dimension = 2;
geom_in.GeometryType = 'cart';
geom_in.xmin_bound_type = glob.IncidentIsotropic;
geom_in.xmax_bound_type = glob.Vacuum;
geom_in.ymin_bound_type = glob.Reflecting;
geom_in.ymax_bound_type = glob.Reflecting;
geom_in.zmin_bound_type = glob.Reflecting;
geom_in.zmax_bound_type = glob.Reflecting;
% ---
dat_in.SpatialMethod = 'PWLD';
dat_in.FEMDegree = 1;
dat_in.FEMLumping = false;
% ---
dat_in.QuadType = 'LS';
dat_in.SnLevels = 4;
dat_in.AzimuthalLevels = 18;
dat_in.PolarLevels = 4;
dat_in.QuadAngles  = [1,1]/norm([1,1]);  % Angles for manual set
dat_in.QuadWeights = 1;                  % Weights for manual set
% ---
dat_in.refineMesh = 0;
dat_in.refinementLevels = 0;
dat_in.AMRIrregularity = 0;
dat_in.refinementTolerance = 0.0;
dat_in.projectSolution = 0;
% ---
dat_in.TotalXS = 1;
dat_in.RHSFunc = {@ZeroTransportFunction};
dat_in.SolFunc = {@ExactSol_IncidentIsotropicLeftFace};
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input(dat_in, geom_in);
% Modify path
if strcmpi(dat_in.GeometryType, 'cart')
    gt = 'Cartesian';
elseif strcmpi(dat_in.GeometryType, 'tri')
    gt = 'Triangular';
elseif strcmpi(dat_in.GeometryType, 'poly')
    gt = 'Polygonal';
end



% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
[data, sol, ~, ~, ~] = execute_problem(data, geometry);

