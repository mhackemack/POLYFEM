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
out.data.TotalXS      = data.Neutronics.TotalXS;
out.data.AbsorbXS     = data.Neutronics.AbsorbXS;
out.data.ScatteringXS = data.Neutronics.ScatteringXS;
out.data.DiffusionXS  = data.Neutronics.DiffusionXS;
out.data.BCFlags      = data.Neutronics.BCFlags;
out.data.BCVals       = data.Neutronics.BCVals;
out.data.IP_Constant  = data.Neutronics.IP_Constant;
out.data.AccelType    = data.Neutronics.AccelType;
out.data.EnergyShape  = data.Neutronics.EnergyShape;
% Average Cross Sections
out.data.numberEnergyGroups = data.Neutronics.numberEnergyGroups;
out.data.AveTotalXS         = data.Neutronics.AveTotalXS;
out.data.AveAbsorbXS        = data.Neutronics.AveAbsorbXS;
out.data.AveScatteringXS    = data.Neutronics.AveScatteringXS;
out.data.AveDiffusionXS     = data.Neutronics.DiffusionXS;
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
% Scattering/Projection Matrix
% ------------------------------------------------------------------------------
out.ScatteringMatrix = build_0th_scattering_matrices(data,out.mesh,out.dof,out.fe);
out.ProjectionMatrix = build_projection_operator(out);
% Quadrature Inputs
% ------------------------------------------------------------------------------
out.Quadrature = inputs.quadrature{q};
% Offset Information
% ------------------------------------------------------------------------------
out.offset = inputs.phase{m}.offset;
% Phase Inputs
% ------------------------------------------------------------------------------
p = inputs.phase{m};