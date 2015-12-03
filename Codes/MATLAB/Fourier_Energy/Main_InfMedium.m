%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          XS Fourier Analysis Script
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
clear; clc; close all; format long e
% Clear Project Space
% ------------------------------------------------------------------------------
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
global glob; glob = get_globals('Home');
inp = '';
addpath(['inputs/InfMedium_AllComponents',inp]); % This one must be last to properly switch input files
% ------------------------------------------------------------------------------
