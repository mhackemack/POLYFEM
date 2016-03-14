%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Spectral Radius Transport Run Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Project Space
% ------------------------------------------------------------------------------
clc; close all; format long e; clear;
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Office');
inp = 'Transport_Quad_Sol';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, ~] = load_user_input();
% Begin user input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BF_names = {'PWLD'};
% BF_names = {'PWLD','WACHSPRESS','MV','MAXENT'};
BF_orders = [2];
data.problem.saveVTKSolution = 1;
data.problem.Dimension = 2;
% End user input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartesian Meshes
% ------------------------------------------------------------------------------
x = linspace(0,1,2);
geometry = CartesianGeometry(2,x,x);
for b=1:length(BF_names)
    now_name = BF_names{b};
    data.Neutronics.SpatialMethod = BF_names{b};
    data.problem.Name = sprintf('cart_%s_k2',now_name);
    [data, geometry] = process_input_data(data, geometry);
    data = cleanup_neutronics_input_data(data, geometry);
    [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
end
