%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Transport Method of Manufactured Solutions (MMS) Script
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
% ------------------------------------------------------------------------------
clc; close all; format long e
fpath = get_path(); 
addpath(fpath);
% Populate global space
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Office'); glob.print_info = false;
inp = 'Transport_MMS';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Define Problem Space to Operate on
% ------------------------------------------------------------------------------
g_type = 'cart';

% Load generic data
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input();
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
[data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
