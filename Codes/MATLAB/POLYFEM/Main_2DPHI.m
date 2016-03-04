%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D PHI Problem Run Script
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
inp = '2DPHI';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
sigt = [10,20,40,80,160,320,640];
c    = [0.9,0.99,0.999,0.9999,0.99999,0.999999];
nstrips = 12;
% Load data and perform error checking
% ------------------------------------------------------------------------------

% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, sol, geometry, DoF, FE] = execute_problem(data, geometry);

