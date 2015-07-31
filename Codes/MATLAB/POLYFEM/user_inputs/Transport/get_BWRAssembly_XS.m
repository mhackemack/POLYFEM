%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          BWR Assembly Benchmark Cross-Sections
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate all cross section data needed for
%                   the BWR Assembly benchmark cases (ANL13).
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_BWRAssembly_XS( data )
% Geometry Information
% --------------------
data.problem.NumberMaterials = 7;
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
data.Neutronics.Transport.TotalXS(1,:) =   [2.531e-1, 5.732e-1];
data.Neutronics.Transport.AbsorbXS(1,:) =  [8.983e-3, 5.892e-2];
data.Neutronics.Transport.FissionXS(1,:) = [2.281e-3, 4.038e-2];
data.Neutronics.Transport.NuBar(1,:) =     [5.925e-3, 9.817e-2]./[2.281e-3, 4.038e-2];
data.Neutronics.Transport.FissSpec(1,:) =  [1, 0];
data.Neutronics.Transport.ScatteringXS(1,1,2,1) = 1.069e-2;
% Material 2
data.Neutronics.Transport.TotalXS(2,:) =   [2.536e-1, 5.767e-1];
data.Neutronics.Transport.AbsorbXS(2,:) =  [8.726e-3, 5.174e-2];
data.Neutronics.Transport.FissionXS(2,:) = [2.003e-3, 3.385e-2];
data.Neutronics.Transport.NuBar(2,:) =     [5.242e-3, 8.228e-2]./[2.003e-3, 3.385e-2];
data.Neutronics.Transport.FissSpec(2,:) =  [1, 0];
data.Neutronics.Transport.ScatteringXS(2,1,2,1) = 1.095e-2;
% Material 3
data.Neutronics.Transport.TotalXS(3,:) =   [2.535e-1, 5.797e-1];
data.Neutronics.Transport.AbsorbXS(3,:) =  [8.587e-3, 4.717e-2];
data.Neutronics.Transport.FissionXS(3,:) = [1.830e-3, 2.962e-2];
data.Neutronics.Transport.NuBar(3,:) =     [4.820e-3, 7.200e-2]./[1.830e-3, 2.962e-2];
data.Neutronics.Transport.FissSpec(3,:) =  [1, 0];
data.Neutronics.Transport.ScatteringXS(3,1,2,1) = 1.112e-2;
% Material 4
data.Neutronics.Transport.TotalXS(4,:) =   [2.533e-1, 5.837e-1];
data.Neutronics.Transport.AbsorbXS(4,:) =  [8.480e-3, 4.140e-2];
data.Neutronics.Transport.FissionXS(4,:) = [1.632e-3, 2.428e-2];
data.Neutronics.Transport.NuBar(4,:) =     [4.337e-3, 5.900e-2]./[1.632e-3, 2.428e-2];
data.Neutronics.Transport.FissSpec(4,:) =  [1, 0];
data.Neutronics.Transport.ScatteringXS(4,1,2,1) = 1.113e-2;
% Material 5
data.Neutronics.Transport.TotalXS(5,:) =   [2.506e-1, 5.853e-1];
data.Neutronics.Transport.AbsorbXS(5,:) =  [9.593e-3, 1.626e-1];
data.Neutronics.Transport.FissionXS(5,:) = [2.155e-3, 9.968e-3];
data.Neutronics.Transport.NuBar(5,:) =     [5.605e-3, 2.424e-2]./[2.155e-3, 9.968e-3];
data.Neutronics.Transport.FissSpec(5,:) =  [1, 0];
data.Neutronics.Transport.ScatteringXS(5,1,2,1) = 1.016e-2;
% Material 6
data.Neutronics.Transport.TotalXS(6,:) =   [2.172e-1, 4.748e-1];
data.Neutronics.Transport.AbsorbXS(6,:) =  [1.043e-3, 4.394e-3];
data.Neutronics.Transport.FissionXS(6,:) = [0, 0];
data.Neutronics.Transport.NuBar(6,:) =     [0, 0];
data.Neutronics.Transport.FissSpec(6,:) =  [0, 0];
data.Neutronics.Transport.ScatteringXS(6,1,2,1) = 9.095e-3;
% Material 7
data.Neutronics.Transport.TotalXS(7,:) =   [2.476e-1, 1.123e+0];
data.Neutronics.Transport.AbsorbXS(7,:) =  [1.983e-4, 7.796e-3];
data.Neutronics.Transport.FissionXS(7,:) = [0, 0];
data.Neutronics.Transport.NuBar(7,:) =     [0, 0];
data.Neutronics.Transport.FissSpec(7,:) =  [0, 0];
data.Neutronics.Transport.ScatteringXS(7,1,2,1) = 3.682e-2;
% Fixup Scattering Cross Sections
data = fixup_sxs(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = fixup_sxs( data )
for m=1:7
    tx = data.Neutronics.Transport.TotalXS(m,:) - data.Neutronics.Transport.AbsorbXS(m,:);
    data.Neutronics.Transport.ScatteringXS(m,1,1,1) = tx(1) - data.Neutronics.Transport.ScatteringXS(m,1,2,1);
    data.Neutronics.Transport.ScatteringXS(m,2,2,1) = tx(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%