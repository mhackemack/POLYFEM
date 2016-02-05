%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Combine Data Structures
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
function [out, p] = combine_input_set(data, inputs, m, q)
% Data Inputs
% ------------------------------------------------------------------------------
out.data.TotalXS      = data.Neutronics.Transport.TotalXS;
out.data.AbsorbXS     = data.Neutronics.Diffusion.AbsorbXS;
out.data.ScatteringXS = data.Neutronics.Transport.ScatteringXS;
out.data.DiffusionXS  = data.Neutronics.Diffusion.DiffusionXS;
out.data.BCFlags      = data.Neutronics.Transport.BCFlags;
out.data.BCVals       = data.Neutronics.Transport.BCVals;
out.data.IP_Constant  = data.Neutronics.IP_Constant;
out.data.AccelType    = data.Neutronics.AccelType;
% Average Cross Sections
out.data.numberEnergyGroups = data.Neutronics.numberEnergyGroups;
out.data.AveTotalXS         = data.Neutronics.Transport.AveTotalXS;
out.data.AveAbsorbXS        = data.Neutronics.Diffusion.AveAbsorbXS;
out.data.AveScatteringXS    = data.Neutronics.Transport.AveScatteringXS;
out.data.AveDiffusionXS     = data.Neutronics.Diffusion.DiffusionXS;
% Hybrid Transport Inputs
% ------------------------------------------------------------------------------
if strcmp(data.Neutronics.Transport.transportType, 'hybrid')
    out.data.StabilizationMethod = data.Neutronics.Transport.StabilizationMethod;
    out.data.FluxStabilization = data.Neutronics.Transport.FluxStabilization;
    out.data.CurrentStabilization = data.Neutronics.Transport.CurrentStabilization;
    if strcmp(out.data.StabilizationMethod, 'upwind')
        out.data.StabilizationType = 0;
    elseif strcmp(out.data.StabilizationMethod, 'EGDG')
        out.data.StabilizationType = 1;
    elseif strcmp(out.data.StabilizationMethod, 'LDG')
        out.data.StabilizationType = 2;
    elseif strcmp(out.data.StabilizationMethod, 'modified')
        out.data.StabilizationType = 3;
    end
end
% Geometry Inputs
% ------------------------------------------------------------------------------
out.mesh = inputs.meshes{m};
out.dof = inputs.dofs{m};
out.fe = inputs.fes{m};
% Scattering Matrix
% ------------------------------------------------------------------------------
out.ScatteringMatrix = build_0th_scattering_matrices(data,out.mesh,out.dof,out.fe);
% Quadrature Inputs
% ------------------------------------------------------------------------------
out.Quadrature = inputs.quadrature{q};
% Offset Information
% ------------------------------------------------------------------------------
out.offset = inputs.phase{m}.offset;
% Phase Inputs
% ------------------------------------------------------------------------------
p = inputs.phase{m};