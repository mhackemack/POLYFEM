%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate MMS Error
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to compute the L2 MMS Error.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err_out = calculate_MMS_error(data, mesh, DoF, FE, sol)
% Get Basic Info and set Output Info
% ----------------------------------
if strcmp(data.Neutronics.transportMethod, 'Diffusion')
    ndat = data.Neutronics.Diffusion;
elseif strcmp(data.Neutronics.transportMethod, 'Transport')
    ndat = data.Neutronics.Transport;
end
if ~iscell(sol), sol = {sol}; end
nout = size(sol,1);
err_out = zeros(nout,1);
% Loop through Solutions
% ----------------------
for i=1:nout
    % Loop through Cells
    for c=1:mesh.TotalCells
        cdofs = DoF.ConnectivityArray{c};
        tsol = sol{i}(cdofs);
        qw = FE.CellQuadWeights{c};
        qx = FE.CellQuadNodes{c};
        cb = FE.CellBasisValues{c};
        sval = cb*tsol;
        fval = ndat.ExactSolution{i}(qx);
        val = qw'*((fval - sval).*(fval-sval));
        err_out(i) = err_out(i) + val;
    end
    err_out(i) = sqrt(err_out(i));
end
clearvars -except err_out;