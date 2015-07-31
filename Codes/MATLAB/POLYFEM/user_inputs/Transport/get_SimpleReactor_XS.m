%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Simple Reactor XS File
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate all cross section data needed for
%                   a simple reactor problem with 2 energy groups and 2
%                   materials - 1 fission material and water.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_SimpleReactor_XS( data )
% Geometry Information
% --------------------
data.problem.NumberMaterials = 2;
% General Neutronics Information
% ------------------------------
data.Neutronics.numberEnergyGroups = 2;
data.Neutronics.Transport.fluxMoments = 0;
% Neutronics Transport Cross-Sections
% -----------------------------------
% Allocate Cross-Sections
nm = data.problem.NumberMaterials;
ng = data.Neutronics.numberEnergyGroups;
nf = data.Neutronics.Transport.fluxMoments + 1;
data.Neutronics.Transport.TotalXS = zeros(nm, ng);
data.Neutronics.Transport.AbsorbXS = zeros(nm, ng);
data.Neutronics.Transport.FissionXS = zeros(nm, ng);
data.Neutronics.Transport.NuBar = zeros(nm, ng);
data.Neutronics.Transport.FissSpec = zeros(nm, ng);
data.Neutronics.Transport.ExtSource = zeros(nm, ng);
data.Neutronics.Transport.ScatteringXS = zeros(nm,ng,ng,nf);
% Material 1
data.Neutronics.Transport.TotalXS(1,:) =   [2.476e-1, 1.123e+0];
data.Neutronics.Transport.AbsorbXS(1,:) =  [1.983e-4, 7.796e-3];
data.Neutronics.Transport.FissionXS(1,:) = [0, 0];
data.Neutronics.Transport.NuBar(1,:) =     [0, 0];
data.Neutronics.Transport.FissSpec(1,:) =  [0, 0];
data.Neutronics.Transport.ScatteringXS(1,1,2,1) = 3.682e-2;
% Material 2
data.Neutronics.Transport.TotalXS(2,:) =   [2.531e-1, 5.732e-1];
data.Neutronics.Transport.AbsorbXS(2,:) =  [8.983e-3, 5.892e-2];
data.Neutronics.Transport.FissionXS(2,:) = [2.281e-3, 4.038e-2];
data.Neutronics.Transport.NuBar(2,:) =     [5.925e-3, 9.817e-2]./[2.281e-3, 4.038e-2];
data.Neutronics.Transport.FissSpec(2,:) =  [1, 0];
data.Neutronics.Transport.ScatteringXS(2,1,2,1) = 1.069e-2;
% Fixup Scattering Cross Sections
data = fixup_sxs(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = fixup_sxs( data )
for m=1:data.problem.NumberMaterials
    tx = data.Neutronics.Transport.TotalXS(m,:) - data.Neutronics.Transport.AbsorbXS(m,:);
    data.Neutronics.Transport.ScatteringXS(m,1,1,1) = tx(1) - data.Neutronics.Transport.ScatteringXS(m,1,2,1);
    data.Neutronics.Transport.ScatteringXS(m,2,2,1) = tx(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%