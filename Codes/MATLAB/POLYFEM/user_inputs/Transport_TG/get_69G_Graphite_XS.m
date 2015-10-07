%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          69 Group Graphite Cross Sections
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate all cross section data needed for
%                   the WIMS 69 group graphite cross section library.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_69G_Graphite_XS( data, ord )
% Set Some Information
% --------------------
data.problem.NumberMaterials = 1;
data.Groups.NumberEnergyGroups = 69;
data.Groups.FastGroups = 1:28;
data.Groups.ThermalGroups = 29:69;
% Neutronics Transport Cross-Sections
% -----------------------------------
% Allocate Cross-Sections
nm = data.problem.NumberMaterials;
ng = data.Groups.NumberEnergyGroups;
nf = ord + 1;
data.XS(1).TotalXS = zeros(nm, ng);
data.XS(1).AbsorbXS = zeros(nm, ng);
data.XS(1).FissionXS = zeros(nm, ng);
data.XS(1).NuBar = zeros(nm, ng);
data.XS(1).FissSpec = zeros(nm, ng);
data.XS(1).ExtSource = zeros(nm, ng);
data.XS(1).ScatteringXS = zeros(nm,ng,ng,nf);
% Get Total Cross Sections
load('user_inputs/Transport_TG/69G_graphite/MT_1.mat');
data.XS(1).TotalXS(1,:) = mat; clear mat;
% Get Scattering Cross Sections
load('user_inputs/Transport_TG/69G_graphite/MT_2500.mat');
data.XS(1).ScatteringXS(1,:,:,1:nf) = mat(:,:,1:nf); clear mat;
% Build Absorption Cross Sections
data.XS(1).AbsorbXS(1,:) = data.XS(1).TotalXS(1,:);
for g=1:ng
    for gg=1:ng
        data.XS(1).AbsorbXS(1,g) = data.XS(1).AbsorbXS(1,g) - data.XS(1).ScatteringXS(1,gg,g,1);
    end
end