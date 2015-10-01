%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set External Diffusion Source
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function src = SetExtSource_Diffusion(data,xsid,groups,mesh,DoF,FE)
% Get some preliminary information
% ------------------------------------------------------------------------------
ngroups = length(groups);
XS = data.XS(xsid);
% Allocate memory space
% ------------------------------------------------------------------------------
src = cell(ngroups,1);
for g=1:ngroups
    src{g} = zeros(DoF.TotalDoFs,1);
end
% Exit with vectors of zeros if there is no inscattering
if ~XS.HasExtSource, return; end
% Build the MMS source
% ------------------------------------------------------------------------------
if data.MMS.PerformMMS
    % Loop through cells
    for c=1:mesh.TotalCells
        cn = DoF.ConnectivityArray{c};
        qx = FE.CellQuadNodes{c};
        qw = FE.CellQuadWeights{c};
        cb = FE.CellBasisValues{c};
        % Loop through energy groups
        for g=1:ngroups
            gfunc = XS.ExtSource{groups(g)};
            src{g}(cn) = src{g}(cn) + cb'*(qw.*gfunc(qx));
        end
    end
end
% Build the constant distributed source
% ------------------------------------------------------------------------------
if ~data.MMS.PerformMMS
    % Loop through cells
    for c=1:mesh.TotalCells
        cmat = mesh.MatID(c);
        cn = DoF.ConnectivityArray{c};
        F = FE.CellFunctionMatrix{c};
        % Loop through energy groups
        for g=1:ngroups
            src{g}(cn) = src{g}(cn) + XS.ExtSource(cmat,groups(g))*F;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%