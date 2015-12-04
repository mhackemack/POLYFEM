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
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Home');
inp = 'InfMedium_AllComponents';
addpath(['inputs/',inp]); % This one must be last to properly switch input files
% ------------------------------------------------------------------------------
data = load_user_input();
% ------------------------------------------------------------------------------
% Loop through materials to analyze and perform analysis
for m=1:data.NumberMaterialsToAnalyze
    tmat = data.Materials{m};
    % Loop through components within the material
    for c=1:data.Materials{m}.NumberComponents
        
    end
end