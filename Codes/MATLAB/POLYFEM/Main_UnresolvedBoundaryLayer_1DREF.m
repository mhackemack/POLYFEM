%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D Unresolved Boundary Layer (grazing angle and opposing
%                   reflecting y-dimension boundaries)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
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
clc; close all; format long e; clear;
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Office');
inp = 'Transport_UBL_1DREF';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
sdm = {'LAGRANGE'};
fedeg = [2,3];
dat_in.FEMLumping = 0;
% ---
dat_in.QuadType = 'PGLC';
dat_in.SnLevels = 16;
dat_in.AzimuthalLevels = 1;
dat_in.PolarLevels = 16;
dat_in.PolarDimension = 1;
% ---
geom_in.Dimension = 1;
geom_in.GeometryType = 'cart';
geom_in.Lx = 1;
geom_in.x = unique([linspace(0,0.01,1001),linspace(0.01,1,1001)]);
% geom_in.x = unique([linspace(0,0.01,201),linspace(0.01,.1,201),linspace(.1,1,101)]);
geom_in.xmin_bound_type = glob.Function;
geom_in.xmax_bound_type = glob.Vacuum;
geom_in.ymin_bound_type = glob.Reflecting;
geom_in.ymax_bound_type = glob.Reflecting;
geom_in.xmin_val = str2func(sprintf('BoundaryFunc_UBL_S%d',dat_in.PolarLevels));
geom_in.xmax_val = 0;
geom_in.ymin_val = 0;
geom_in.ymax_val = 0;
% ---
dat_in.refineMesh = 0;
dat_in.refinementLevels = 4;
dat_in.AMRIrregularity = 1;
dat_in.refinementTolerance = 0.3;
dat_in.projectSolution = 1;
% ---
sigt = [500];
dat_in.RHSFunc = {@ZeroTransportFunction};
% Execute Problem Suite
% ------------------------------------------------------------------------------
% Loop through FEM degrees
for k=1:length(fedeg)
    dat_in.FEMDegree = fedeg(k);
    % Loop through basis functions
    for s=1:length(sdm)
        dat_in.SpatialMethod = sdm{s};
        % Loop through total cross sections
        for t=1:length(sigt)
            dat_in.TotalXS = sigt(t);
            data = load_user_input(dat_in, geom_in);
            [data,geometry] = load_geometry_input(data, geom_in);
            data.problem.Path = 'Transport/UnresolvedBoundaryLayer_1DREF';
            if dat_in.FEMLumping
                data.problem.Name = sprintf('S%d_Lx%d_sig%d_L%s%d',dat_in.SnLevels,geom_in.Lx,sigt(t),sdm{s},fedeg(k));
            else
                data.problem.Name = sprintf('S%d_Lx%d_sig%d_U%s%d',dat_in.SnLevels,geom_in.Lx,sigt(t),sdm{s},fedeg(k));
            end
            % Run problem iteration
            [data, geometry] = process_input_data(data, geometry);
            data = cleanup_neutronics_input_data(data, geometry);
            [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
            plot_solution(geometry,DoF,FE,sol.flux);
        end
    end
end

