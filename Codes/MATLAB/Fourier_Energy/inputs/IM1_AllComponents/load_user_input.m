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
data.OutputName = 'IM1_Materials';
data.Energy.NumberEngeryGroups = 99;
data.Energy.PnOrder = 8;
data.Energy.FastGroups = 1:42;
data.Energy.ThermalGroups = 43:99;
% ------------------------------------------------------------------------------
p = 1;
% data.Materials{p}.MaterialName = 'Al_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Al_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Am241_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Am241_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Ar40_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Ar40_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'B10_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'B10_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'B11_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'B11_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Be9_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Be9_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Cr52_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Cr52_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'F19_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'F19_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Fe56_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Fe56_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'FG_CNat_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'FG_CNat_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'FG_H1_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'FG_H1_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'graphite_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'graphite_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'H2O_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'H2O_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Mn55_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Mn55_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'N14_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'N14_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Ni58_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Ni58_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'Ni60_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'Ni60_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'O16_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'O16_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
% data.Materials{p}.MaterialName = 'PolyH1_99G';
% data.Materials{p}.NumberComponents = 1;
% data.Materials{p}.ComponentNames = {'PolyH1_99G'};
% data.Materials{p}.ComponentDensities = 1.0;
% % ------------------------------------------------------------------------------
% p = p + 1;
data.Materials{p}.MaterialName = 'IM1_Air';
data.Materials{p}.NumberComponents = 4;
data.Materials{p}.ComponentNames = {'FG_CNat_99G','N14_99G','O16_99G','Ar40_99G'};
data.Materials{p}.ComponentDensities = [7.4906E-9,3.9123E-5,1.0511E-5,2.3297E-7];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_Wood';
data.Materials{p}.NumberComponents = 3;
data.Materials{p}.ComponentNames = {'FG_H1_99G','FG_CNat_99G','O16_99G'};
data.Materials{p}.ComponentDensities = [2.0752E-2,1.4520E-2,1.0376E-2];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_Boral';
data.Materials{p}.NumberComponents = 4;
data.Materials{p}.ComponentNames = {'Al_99G','B10_99G','B11_99G','FG_CNat_99G'};
data.Materials{p}.ComponentDensities = [3.8193E-2,7.1036E-3,2.8593E-2,8.9241E-3];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_BHDPE';
data.Materials{p}.NumberComponents = 4;
data.Materials{p}.ComponentNames = {'PolyH1_99G','FG_CNat_99G','B10_99G','B11_99G'};
data.Materials{p}.ComponentDensities = [5.0859E-2,2.5429E-2,6.6256E-3,2.6669E-2];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_HDPE';
data.Materials{p}.NumberComponents = 2;
data.Materials{p}.ComponentNames = {'PolyH1_99G','FG_CNat_99G'};
data.Materials{p}.ComponentDensities = [8.1570E-2,4.0787E-2];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_AmBe';
data.Materials{p}.NumberComponents = 3;
data.Materials{p}.ComponentNames = {'Am241_99G','Be9_99G','O16_99G'};
data.Materials{p}.ComponentDensities = [1.1649E-3,1.9077E-1,2.3298E-3];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_Graphite';
data.Materials{p}.NumberComponents = 2;
data.Materials{p}.ComponentNames = {'graphite_99G','B10_99G'};
data.Materials{p}.ComponentDensities = [8.5238E-2,2.4335449e-06];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_Steel';
data.Materials{p}.NumberComponents = 4;
data.Materials{p}.ComponentNames = {'Cr52_99G','Mn55_99G','Fe56_99G','Ni58_99G'};
data.Materials{p}.ComponentDensities = [1.7428E-2,1.7363E-3,5.9358E-2,7.7199E-3];
% ------------------------------------------------------------------------------
p = p + 1;
data.Materials{p}.MaterialName = 'IM1_BF3';
data.Materials{p}.NumberComponents = 3;
data.Materials{p}.ComponentNames = {'B10_99G','B11_99G','F19_99G'};
data.Materials{p}.ComponentDensities = [6.4458E-6,2.6858E-7,2.0143E-5];
% ------------------------------------------------------------------------------
data.Materials = data.Materials';
data.NumberMaterialsToAnalyze = length(data.Materials);
