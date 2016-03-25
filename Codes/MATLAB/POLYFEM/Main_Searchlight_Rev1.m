%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Searchlight Problem Run Script (Rev1)
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
inp = 'Searchlight_Aligned';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input();
% Modify path
data.problem.Path = sprintf('Transport/Searchlight_Aligned/%s_k%d',data.Neutronics.SpatialMethod,data.Neutronics.FEMDegree);
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, ~, ~, ~, ~] = execute_problem(data, geometry);
% Process Outputs
% ------------------------------------------------------------------------------
adir = data.Neutronics.Transport.QuadAngles';
% Build data storage structures
nr = data.problem.refinementLevels + 1;
dofnum = zeros(nr,1);
