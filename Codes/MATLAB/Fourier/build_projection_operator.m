%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Preconditioner Projection Operation
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = build_projection_operator(input)
% Retrieve some preliminary information
% ------------------------------------------------------------------------------
ng = input.data.numberEnergyGroups; ndofs = input.dof.TotalDoFs;
ntot = ndofs*ng;
eshape = input.data.ErrorShape;
% Allocate memory space
% ------------------------------------------------------------------------------
g_offset = (1:ng)*ndofs - ndofs;
% P = zeros(ndofs,ndofs,ng);
P = zeros(ntot,ndofs);
% Loop through mesh cells and build projection operator
% ------------------------------------------------------------------------------
for c=1:input.mesh.TotalCells
    mat = input.mesh.MatID(c);
    cn  = input.dof.ConnectivityArray{c}; ncn = length(cn);
    tmat = diag(ones(ncn,1));
    for g=1:ng
%         P(cn,cn,g) = eshape(mat,g)*tmat;
        gcn = cn + g_offset(g);
        P(gcn,cn) = eshape(mat,g)*tmat;
    end
end