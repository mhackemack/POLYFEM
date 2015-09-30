%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set Total Transport Source
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function src = SetSource_Transport(XS, mquad, groups, flux, mesh, DoF, FE)
% Get some preliminary information
% ------------------------------------------------------------------------------
ngrps = length(groups);
% Allocate memory space
% ------------------------------------------------------------------------------

% Loop through spatial cells and build source
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    
    % Double loop through energy groups
    
end