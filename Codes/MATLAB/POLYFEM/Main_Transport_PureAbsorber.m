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
dat_in.Name = 'LeftIncidentIsotropic';
% ---
geom_in.Dimension = 2;
geom_in.GeometryType = 'cart';
geom_in.PolyNum = 256;
geom_in.Lx = 1; geom_in.ncellx = 4;
geom_in.Ly = 1; geom_in.ncelly = 4;
geom_in.Lz = 1; geom_in.ncellz = 4;
geom_in.xmin_bound_type = glob.IncidentIsotropic;
geom_in.xmax_bound_type = glob.Vacuum;
geom_in.ymin_bound_type = glob.Reflecting;
geom_in.ymax_bound_type = glob.Reflecting;
geom_in.zmin_bound_type = glob.Reflecting;
geom_in.zmax_bound_type = glob.Reflecting;
geom_in.xmin_val = 2;
geom_in.xmax_val = 0;
geom_in.ymin_val = 0;
geom_in.ymax_val = 0;
geom_in.zmin_val = 0;
geom_in.zmax_val = 0;
% ---
dat_in.SpatialMethod = 'PWLD';
dat_in.FEMDegree = 1;
dat_in.FEMLumping = false;
% ---
dat_in.QuadType = 'LS';
dat_in.SnLevels = 16;
dat_in.AzimuthalLevels = 18;
dat_in.PolarLevels = 4;
dat_in.QuadAngles  = [1,-1]/norm([1,1]);  % Angles for manual set
dat_in.QuadWeights = 1;                  % Weights for manual set
% ---
dat_in.refineMesh = 1;
dat_in.refinementLevels = 3;
dat_in.AMRIrregularity = 1;
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
if strcmpi(geom_in.GeometryType, 'cart')
    gt = 'Cartesian';
elseif strcmpi(geom_in.GeometryType, 'tri')
    gt = 'Triangular';
elseif strcmpi(geom_in.GeometryType, 'poly')
    gt = 'Polygonal';
end
data.problem.Path = ['Transport/PureAbsorber/',gt,'/',dat_in.SpatialMethod,num2str(dat_in.FEMDegree)];
data.problem.Name = dat_in.Name;
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
[data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
% Postprocess Solution Data
% ------------------------------------------------------------------------------
if ~iscell(sol), sol = {sol}; end
% Loop through refinements
for r=1:length(sol)
    
end
