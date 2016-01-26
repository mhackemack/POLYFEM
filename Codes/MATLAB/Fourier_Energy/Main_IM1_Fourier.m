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
global glob; glob = get_globals('Office');
inp = 'IM1_AllComponents';
addpath(['inputs/',inp]); % This one must be last to properly switch input files
% ------------------------------------------------------------------------------
data = load_user_input();
% ------------------------------------------------------------------------------
ny = 1e3; x = linspace(1e-8,4*pi,ny);
nm = data.NumberMaterialsToAnalyze;
ng = data.Energy.NumberEngeryGroups;
nmom = data.Energy.PnOrder;
y_P0_noaccel = zeros(ny,nm);
y_P0_accel   = zeros(ny,nm);
y_P1_noaccel = zeros(ny,nm);
y_P1_accel   = zeros(ny,nm);
% ------------------------------------------------------------------------------
func_P0_noaccel = get_2G_fourier_func('unaccelerated',0);
func_P0_accel   = get_2G_fourier_func('accelerated',0);
func_P1_noaccel = get_2G_fourier_func('unaccelerated',1);
func_P1_accel   = get_2G_fourier_func('accelerated',1);
% ------------------------------------------------------------------------------
mat_str_out = cell(1,data.NumberMaterialsToAnalyze);
% ------------------------------------------------------------------------------
% Loop through materials to analyze and perform analysis
for m=1:data.NumberMaterialsToAnalyze
    tmat = data.Materials{m};
    matname = data.Materials{m}.MaterialName;
    mat_str_out{m} = matname;
    tcompnames = data.Materials{m}.ComponentNames;
    tcompdens = data.Materials{m}.ComponentDensities;
    % Zero out total and scattering cross sections
    T = zeros(ng,ng);
    S = zeros(ng,ng,nmom);
    % Loop through components within the material
    for c=1:data.Materials{m}.NumberComponents
        % Add total xs contribution
        load([glob.XS_path,tcompnames{1},'/MT_1.mat']);
        T = T + tcompdens(c)*diag(mat);
        % Add P0 scattering xs contribution
        load([glob.XS_path,tcompnames{1},'/MT_2500.mat']);
        S = S + tcompdens(c)*mat;
    end
end