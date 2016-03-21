%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Quadratic Solution Correction
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
function [data, sol] = Quadratic_Solution_Correction(data, mesh, DoF, FE, sol)
% Set the correct solution parameters
% ------------------------------------------------------------------------------
a = 1;
b = 1;
c = 1;
d = 1;
e = 1;
f = 1;
% Retrieve some preliminary information
% ------------------------------------------------------------------------------
ng = data.Neutronics.numberEnergyGroups;
err_out = zeros(ng,1);
flux  = sol.flux;
cflux = flux;
% Loop through mesh cells
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cdofs = DoF.ConnectivityArray{c};
    qw = FE.CellQuadWeights{c};
    qx = FE.CellQuadNodes{c};
    cb = FE.CellBasisValues{c};
    % Loop through energy groups
    for g=1:ng
        tsol = sol{g}(cdofs);
        
    end
end
% Update solutions
% ------------------------------------------------------------------------------
sol.flux = cflux;