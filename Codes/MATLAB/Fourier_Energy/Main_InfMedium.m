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
GSPI_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
GS_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
for m=1:data.NumberMaterialsToAnalyze
    tmat = data.Materials{m};
    matname = data.Materials{m}.MaterialName;
    GSPI_str_out{1,m} = matname;
    GS_str_out{1,m} = matname;
    tcompnames = data.Materials{m}.ComponentNames;
    tcompdens = data.Materials{m}.ComponentDensities;
    % Get energy bounds from first component in material
    load([glob.XS_path,tcompnames{1},'/Energy_Bounds.mat']); E = mat;
    Eave = (E(1:99) + E(2:100))./2; Ediff = E(1:99) - E(2:100);
    % Zero out total and scattering cross sections
    T = zeros(data.Energy.NumberEngeryGroups);
    S0 = zeros(data.Energy.NumberEngeryGroups);
    % Loop through components within the material
    for c=1:data.Materials{m}.NumberComponents
        % Add total xs contribution
        load([glob.XS_path,tcompnames{1},'/MT_1.mat']);
        T = T + tcompdens(c)*diag(mat);
        % Add P0 scattering xs contribution
        load([glob.XS_path,tcompnames{1},'/MT_2500.mat']);
        S0 = S0 + tcompdens(c)*mat(:,:,1);
    end
    % Generate Infinite Medium EigenSpectrums
    display(sprintf('Computing Spectrum for material %d of %d.',m,data.NumberMaterialsToAnalyze))
    tg = data.Energy.ThermalGroups;
    [PeJ,PeGS,eJ,eGS] = calculate_inf_medium_P0_eigenspectrum(tg,T,S0);
    PeGS.EigenSpectrum = PeGS.EigenSpectrum / sum(PeGS.EigenSpectrum);
    [~,mi] = max(abs(eGS.EigenValue));
    eGS_ES = eGS.EigenSpectrum(:,mi) / sum(eGS.EigenSpectrum(:,mi));
    % Loop through thermal groups for output
    for g=1:length(data.Energy.ThermalGroups)
        GSPI_str_out{g+1,m} = num2str(PeGS.EigenSpectrum(g),'%16.8e');
        GS_str_out{g+1,m} = num2str(eGS_ES(g),'%16.8e');
    end
end
cell2csv('GSPI_Out',GSPI_str_out);
cell2csv('GS_Out',GS_str_out);
% ------------------------------------------------------------------------------