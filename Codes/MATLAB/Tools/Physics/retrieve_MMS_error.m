%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Retrieve MMS Error
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to retrieve MMS error values from data
%                   structures.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dofs, err] = retrieve_MMS_error(data, sol)
ref_bool = data.problem.refineMesh;
ng = data.Neutronics.numberEnergyGroups;
% Allocate Memory
if ~ref_bool
    dofs = sol.SpatialDoFs;
    err = sol.MMS_error;
    return
end
ref_lvls = data.problem.refinementLevels;
dofs = zeros(ref_lvls+1,1);
err = zeros(ref_lvls+1,ng);
% Retrieve dofs and error
for i=1:ref_lvls+1
    dofs(i) = sol{i}.SpatialDoFs;
    err(i,:) = sol{i}.MMS_error;
end