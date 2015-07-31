%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          IAEA-EIR-2 Benchmark Cross Sections
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate all cross section data needed for
%                   the IAEA-EIR-2 benchmark cases.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_Yaqi_XS( data )
% Geometry Information
% --------------------
data.problem.NumberMaterials = 4;
% General Neutronics Information
% ------------------------------
data.Neutronics.numberEnergyGroups = 1;
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
data.Neutronics.Transport.TotalXS(1,:) =   1.00;
data.Neutronics.Transport.AbsorbXS(1,:) =  1.00;
data.Neutronics.Transport.FissionXS(1,:) = 0.00;
data.Neutronics.Transport.NuBar(1,:) =     0.00;
data.Neutronics.Transport.FissSpec(1,:) =  0.00;
data.Neutronics.Transport.ExtSource(1,:) = 1.00;
data.Neutronics.Transport.ScatteringXS(1,:,:,1) = 0.00;
% Material 2
data.Neutronics.Transport.TotalXS(2,:) =   0.001;
data.Neutronics.Transport.AbsorbXS(2,:) =  0.001;
data.Neutronics.Transport.FissionXS(2,:) = 0.00;
data.Neutronics.Transport.NuBar(2,:) =     0.00;
data.Neutronics.Transport.FissSpec(2,:) =  0.00;
data.Neutronics.Transport.ExtSource(2,:) = 0.00;
data.Neutronics.Transport.ScatteringXS(2,:,:,1) = 0.00;
% Material 3
data.Neutronics.Transport.TotalXS(3,:) =   100.00;
data.Neutronics.Transport.AbsorbXS(3,:) =  100.00;
data.Neutronics.Transport.FissionXS(3,:) = 0.00;
data.Neutronics.Transport.NuBar(3,:) =     0.00;
data.Neutronics.Transport.FissSpec(3,:) =  0.00;
data.Neutronics.Transport.ExtSource(3,:) = 0.00;
data.Neutronics.Transport.ScatteringXS(3,:,:,1) = 0.00;
% Material 4
data.Neutronics.Transport.TotalXS(4,:) =   10.0;
data.Neutronics.Transport.AbsorbXS(4,:) =  0.001;
data.Neutronics.Transport.FissionXS(4,:) = 0.00;
data.Neutronics.Transport.NuBar(4,:) =     0.00;
data.Neutronics.Transport.FissSpec(4,:) =  0.00;
data.Neutronics.Transport.ExtSource(4,:) = 0.00;
data.Neutronics.Transport.ScatteringXS(4,:,:,1) = 9.999;