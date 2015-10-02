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
function data = get_Reed_XS( data, xsid )
% Geometry Information
% --------------------
data.problem.NumberMaterials = 5;
% General Neutronics Information
% ------------------------------
data.Groups.NumberEnergyGroups = 1;
data.Transport.PnOrder = 0;
% Neutronics Transport Cross-Sections
% -----------------------------------
% Allocate Cross-Sections
nm = data.problem.NumberMaterials;
ng = data.Groups.NumberEnergyGroups;
nf = data.Transport.PnOrder + 1;
data.XS(xsid).TotalXS = zeros(nm, ng);
data.XS(xsid).AbsorbXS = zeros(nm, ng);
data.XS(xsid).FissionXS = zeros(nm, ng);
data.XS(xsid).NuBar = zeros(nm, ng);
data.XS(xsid).FissSpec = zeros(nm, ng);
data.XS(xsid).ExtSource = zeros(nm, ng);
data.XS(xsid).ScatteringXS = zeros(nm,ng,ng,nf);
% Material 1
data.XS(xsid).TotalXS(1,:) =   50.0;
data.XS(xsid).AbsorbXS(1,:) =  50.0;
data.XS(xsid).FissionXS(1,:) = 0.00;
data.XS(xsid).NuBar(1,:) =     0.00;
data.XS(xsid).FissSpec(1,:) =  0.00;
data.XS(xsid).ExtSource(1,:) = 50.0;
data.XS(xsid).ScatteringXS(1,:,:,1) = 0.00;
% Material 2
data.XS(xsid).TotalXS(2,:) =   5.00;
data.XS(xsid).AbsorbXS(2,:) =  5.00;
data.XS(xsid).FissionXS(2,:) = 0.00;
data.XS(xsid).NuBar(2,:) =     0.00;
data.XS(xsid).FissSpec(2,:) =  0.00;
data.XS(xsid).ExtSource(2,:) = 0.00;
data.XS(xsid).ScatteringXS(2,:,:,1) = 0.00;
% Material 3
data.XS(xsid).TotalXS(3,:) =   0.0001;
data.XS(xsid).AbsorbXS(3,:) =  0.0001;
data.XS(xsid).FissionXS(3,:) = 0.00;
data.XS(xsid).NuBar(3,:) =     0.00;
data.XS(xsid).FissSpec(3,:) =  0.00;
data.XS(xsid).ExtSource(3,:) = 0.00;
data.XS(xsid).ScatteringXS(3,:,:,1) = 0.00;
% Material 4
data.XS(xsid).TotalXS(4,:) =   20.0;
data.XS(xsid).AbsorbXS(4,:) =  0.01;
data.XS(xsid).FissionXS(4,:) = 0.00;
data.XS(xsid).NuBar(4,:) =     0.00;
data.XS(xsid).FissSpec(4,:) =  0.00;
data.XS(xsid).ExtSource(4,:) = 0.10;
data.XS(xsid).ScatteringXS(4,:,:,1) = 19.99;
% Material 5
data.XS(xsid).TotalXS(5,:) =   20.0;
data.XS(xsid).AbsorbXS(5,:) =  0.01;
data.XS(xsid).FissionXS(5,:) = 0.00;
data.XS(xsid).NuBar(5,:) =     0.00;
data.XS(xsid).FissSpec(5,:) =  0.00;
data.XS(xsid).ExtSource(5,:) = 0.00;
data.XS(xsid).ScatteringXS(5,:,:,1) = 19.99;
% Set Diffusion XS
data.XS(xsid).DiffXS =    (1./data.XS(xsid).TotalXS)/3;