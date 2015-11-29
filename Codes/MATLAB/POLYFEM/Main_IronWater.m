%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          IronWater Problem Run Script
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
global glob; glob = get_globals('Home');
inp = 'IronWater';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
dat_in.GeometryType = 'cart';
dat_in.SpatialMethod = 'PWLD';
dat_in.FEMDegree = 1;
% ---
dat_in.QuadType = 'LS';
dat_in.SnLevels = 4;
dat_in.AzimuthalLevels = 14;
dat_in.PolarLevels = 2;
% ---
dat_in.refinementLevels = 26;
dat_in.AMRIrregularity = 1;
dat_in.refinementTolerance = 1/5;
dat_in.projectSolution = 0;
% ---
dat_in.DSASolveMethod = 'PCG';
dat_in.DSAPreconditioner = 'gs';
dat_in.DSATolerance = 1e-3;
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input(dat_in);
% Modify path
if dat_in.projectSolution == 1
    bs = 'WITHBootstrapping';
elseif dat_in.projectSolution == 0
    bs = 'NOBootstrapping';
end
if strcmpi(dat_in.GeometryType, 'cart')
    gt = 'Cartesian';
elseif strcmpi(dat_in.GeometryType, 'tri')
    gt = 'Triangular';
end
data.problem.Path = ['Transport/IronWater/',bs,'/',gt,'/',dat_in.SpatialMethod,num2str(dat_in.FEMDegree)];
% Modify name
if strcmpi(dat_in.QuadType, 'LS')
    oname = [dat_in.QuadType,num2str(dat_in.SnLevels)];
elseif strcmpi(dat_in.QuadType, 'PGLC')
    oname = [dat_in.QuadType,num2str(dat_in.AzimuthalLevels),'-',num2str(dat_in.PolarLevels)];
end
if abs(dat_in.refinementTolerance) < 1e-8
    oname = [oname, '_uniform'];
else
    oname = [oname,sprintf('_Irr=%d_tol=%1.3f',dat_in.AMRIrregularity,dat_in.refinementTolerance)];
end
oname = [oname,'_',dat_in.DSASolveMethod,'_',dat_in.DSAPreconditioner];
data.problem.Name = oname;
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
[data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
% [~, ~, ~, ~, ~] = execute_problem(data, geometry);