%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set PHI Cross Sections
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
function data = set_phi_xs(data, sigt, c)
% Set Cross Sections
sigt = [1/sigt;1];
data.Neutronics.TotalXS = sigt;
data.Neutronics.DiffusionXS = (1/3)./sigt;
data.Neutronics.ScatteringXS = c*sigt;
data.Neutronics.AbsorbXS = (1-c)*sigt;
% Set Average Cross Sections
data.Neutronics.AveTotalXS      = data.Neutronics.TotalXS;
data.Neutronics.AveDiffusionXS  = data.Neutronics.DiffusionXS;
data.Neutronics.AveScatteringXS = data.Neutronics.ScatteringXS;
data.Neutronics.AveAbsorbXS     = data.Neutronics.AbsorbXS;