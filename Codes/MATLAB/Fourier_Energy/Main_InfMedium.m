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
inp = 'InfMedium_GraphiteONLY';
addpath(['inputs/',inp]); % This one must be last to properly switch input files
% ------------------------------------------------------------------------------
data = load_user_input();
ng = data.Energy.NumberEngeryGroups;
% ------------------------------------------------------------------------------
% Loop through materials to analyze and perform analysis
GSPI_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
GS_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
JPI_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
J_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
MJIAPI_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
MJIA_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
MTGPI_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
MTG_str_out = cell(length(data.Energy.ThermalGroups)+1,data.NumberMaterialsToAnalyze);
for m=1:data.NumberMaterialsToAnalyze
    tmat = data.Materials{m};
    matname = data.Materials{m}.MaterialName;
    GSPI_str_out{1,m} = matname; GS_str_out{1,m} = matname;
    JPI_str_out{1,m} = matname; J_str_out{1,m} = matname;
    MTGPI_str_out{1,m} = matname; MTG_str_out{1,m} = matname;
    MJIAPI_str_out{1,m} = matname; MJIA_str_out{1,m} = matname;
    tcompnames = data.Materials{m}.ComponentNames;
    tcompdens = data.Materials{m}.ComponentDensities;
    % Get energy bounds from first component in material
    load([glob.XS_path,tcompnames{1},'/Energy_Bounds.mat']); E = mat;
    Eave = (E(1:ng) + E(2:(ng+1)))./2; Ediff = E(1:ng) - E(2:(ng+1));
    % Zero out total and scattering cross sections
    T = zeros(data.Energy.NumberEngeryGroups);
    S0 = zeros(data.Energy.NumberEngeryGroups);
    % Loop through components within the material
    for c=1:data.Materials{m}.NumberComponents
        % Add total xs contribution
        load([glob.XS_path,tcompnames{c},'/MT_1.mat']);
        T = T + tcompdens(c)*diag(mat);
        % Add P0 scattering xs contribution
        load([glob.XS_path,tcompnames{c},'/MT_2500.mat']);
        S0 = S0 + tcompdens(c)*mat(:,:,1);
    end
    % Generate Infinite Medium EigenSpectrums
    display(sprintf('Computing Spectrum for material %d of %d.',m,data.NumberMaterialsToAnalyze))
    tg = data.Energy.ThermalGroups;
    [PeJ,PeMJIA,PeGS,PeMTG,eJ,eMJIA,eGS,eMTG] = calculate_inf_medium_P0_eigenspectrum(tg,T,S0);
    % Within-group jacobi spectrums
    PeJ.EigenSpectrum = PeJ.EigenSpectrum / sum(PeJ.EigenSpectrum);
    [~,mi] = max(abs(eJ.EigenValue));
    eJ_ES = eJ.EigenSpectrum(:,mi) / sum(eJ.EigenSpectrum(:,mi));
    for g=1:length(data.Energy.ThermalGroups)
        JPI_str_out{g+1,m} = num2str(PeJ.EigenSpectrum(g),'%16.8e');
        J_str_out{g+1,m} = num2str(eJ_ES(g),'%16.8e');
    end
    % MJIA spectrums
    PeMJIA.EigenSpectrum = PeMJIA.EigenSpectrum / sum(PeMJIA.EigenSpectrum);
    [~,mi] = max(abs(eMJIA.EigenValue));
    eMJIA_ES = eMJIA.EigenSpectrum(:,mi) / sum(eMJIA.EigenSpectrum(:,mi));
    for g=1:length(data.Energy.ThermalGroups)
        MJIAPI_str_out{g+1,m} = num2str(PeMJIA.EigenSpectrum(g),'%16.8e');
        MJIA_str_out{g+1,m} = num2str(eMJIA_ES(g),'%16.8e');
    end
    % Gauss-Seidel spectrums
    PeGS.EigenSpectrum = PeGS.EigenSpectrum / sum(PeGS.EigenSpectrum);
    [~,mi] = max(abs(eGS.EigenValue));
    eGS_ES = eGS.EigenSpectrum(:,mi) / sum(eGS.EigenSpectrum(:,mi));
    for g=1:length(data.Energy.ThermalGroups)
        GSPI_str_out{g+1,m} = num2str(PeGS.EigenSpectrum(g),'%16.8e');
        GS_str_out{g+1,m} = num2str(eGS_ES(g),'%16.8e');
    end
    % Modified two-grid spectrums
    PeMTG.EigenSpectrum = PeMTG.EigenSpectrum / sum(PeMTG.EigenSpectrum);
    [~,mi] = max(abs(eMTG.EigenValue));
    eMTG_ES = eMTG.EigenSpectrum(:,mi) / sum(eMTG.EigenSpectrum(:,mi));
    for g=1:length(data.Energy.ThermalGroups)
        MTGPI_str_out{g+1,m} = num2str(PeMTG.EigenSpectrum(g),'%16.8e');
        MTG_str_out{g+1,m} = num2str(eMTG_ES(g),'%16.8e');
    end
end
% ------------------------------------------------------------------------------
% Check if output directory exists
if ~isequal(exist([glob.output_path,data.OutputName], 'dir'),7),mkdir([glob.output_path,data.OutputName]); end
% Gauss-Seidel outputs
% cell2csv([glob.output_path,data.OutputName,'/GSPI_Out'],GSPI_str_out); cell2csv([glob.output_path,data.OutputName,'/GS_Out'],GS_str_out);
% Jacobi outputs
% cell2csv([glob.output_path,data.OutputName,'/JPI_Out'],JPI_str_out); cell2csv([glob.output_path,data.OutputName,'/J_Out'],J_str_out);
% Modified two-grid outputs
% cell2csv([glob.output_path,data.OutputName,'/MTGPI_Out'],MTGPI_str_out); cell2csv([glob.output_path,data.OutputName,'/MTG_Out'],MTG_str_out);
% MJIA outputs
cell2csv([glob.output_path,data.OutputName,'/MJIAPI_Out'],MJIAPI_str_out); cell2csv([glob.output_path,data.OutputName,'/MJIA_Out'],MJIA_str_out);
% ------------------------------------------------------------------------------
