%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          IAEA-EIR-2 Problem Run Script
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
inp = 'EIR2';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
dat_in.GeometryType = 'cart';
dat_in.SpatialMethod = 'PWLD';
dat_in.FEMDegree = 1;
% ---
dat_in.QuadType = 'LS';
dat_in.SnLevels = 8;
dat_in.AzimuthalLevels = 14;
dat_in.PolarLevels = 2;
% ---
dat_in.refinementLevels = 2;
dat_in.AMRIrregularity = 1;
dat_in.refinementTolerance = 0;
dat_in.projectSolution = 0;
% ---
dat_in.DSASolveMethod = 'direct';
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input(dat_in);
% Modify path
if strcmpi(dat_in.GeometryType, 'cart')
    gt = 'Cartesian';
elseif strcmpi(dat_in.GeometryType, 'tri')
    gt = 'Triangular';
end
data.problem.Path = ['Transport/EIR2/',gt,'/',dat_in.SpatialMethod,num2str(dat_in.FEMDegree)];
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
data.problem.Name = oname;
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
[data, sol, ~, ~, ~] = execute_problem(data, geometry);
% Postprocess Solution Data
% ------------------------------------------------------------------------------
ave_flux = [1.5263e0;1.1960e1;5.3968e-1;1.9202e1;8.3364e-1];
mat_map = [5;1;2;3;4];
