%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D/2D/3D Poly Run Script
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
clc; close all; format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Populate global space
% ---------------------
global glob
glob = get_globals('Office');
% Specify User-Specific Input Folder Location
% -------------------------------------------
% inp = 'Diffusion';
inp = 'Transport_Rev1';
% Populate path with additional folders
% -------------------------------------
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Load data and perform error checking
% ------------------------------------
print_heading(now,date);
[data, geometry] = load_user_input();
[data, geometry] = process_input(data, geometry);
% [data, geometry] = process_input_data(data, geometry);
% data = cleanup_neutronics_input_data(data, geometry);
% Execute Problem Suite
% ---------------------
[data, geometry, DoF, FE] = execute_problem_Rev1(data, geometry);
% [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
