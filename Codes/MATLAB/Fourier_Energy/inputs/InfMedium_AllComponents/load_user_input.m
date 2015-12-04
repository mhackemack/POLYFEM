%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          User Input - All Component Infinite Medium
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
function data = load_user_input()
% ------------------------------------------------------------------------------
p = 1;
data.Materials{p}.MaterialName = 'Al_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Al_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Am241_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Am241_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Ar40_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Ar40_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'B10_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'B10_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'B11_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'B11_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Be9_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Be9_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Cr52_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Cr52_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'F19_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'F19_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Fe56_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Fe56_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'FG_CNat_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'FG_CNat_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'FG_H1_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'FG_H1_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'graphite_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'graphite_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'H2O_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'H2O_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Mn55_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Mn55_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'N14_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'N14_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Ni58_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Ni58_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'Ni60_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'Ni60_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'O16_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'O16_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'PolyH1_99G';
data.Materials{p}.NumberComponents = 1;
data.Materials{p}.ComponentNames = {'PolyH1_99G'};
data.Materials{p}.ComponentDensities = 1.0;
% ------------------------------------------------------------------------------
data.Materials = data.Materials';
data.NumberMaterialsToAnalyze = length(data.Materials);
