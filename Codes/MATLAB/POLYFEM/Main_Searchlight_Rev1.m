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
inp = 'Searchlight_Rev1';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
geom_in.Dimension = 2;
geom_in.GeometryType = 'cart';
geom_in.Lx = 1; geom_in.ncellx = 5;
geom_in.Ly = 1; geom_in.ncelly = 5;
% ---
fedeg = [2];
sdm = {'PWLD','MV'};
% ---
dat_in.lvls = 30;
dat_in.irr = 1;
dat_in.tol = 0.1;