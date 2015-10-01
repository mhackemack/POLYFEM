%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set Scattering Transport Source
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function src = SetScatteringSource_Transport(XS,mquad,groups,gin,flux,mesh,DoF,FE)
% Get some preliminary information
% ------------------------------------------------------------------------------
ngroups = length(groups); ngin = length(gin);
% Allocate memory space
% ------------------------------------------------------------------------------
src = cell(ngroups,mquad.TotalFluxMoments);
for g=1:ngroups
    for m=1:mquad.TotalFluxMoments
        src{g,m} = zeros(DoF.TotalDoFs,1);
    end
end
% Exit with vectors of zeros if there is no inscattering
if ngin == 0 || ~XS.HasScattering, return; end
% Loop through spatial cells and build source
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    cn   = DoF.ConnectivityArray{c}; ncn = length(cn);
    M    = FE.CellMassMatrix{c};
    for m=1:mquad.TotalFluxMoments
        k = mquad.MomentOrders(m,1) + 1;
        % Double loop through energy groups
        for g=1:ngroups
            grp = groups(g);
            gvec = zeros(ncn,1);
            for gg=1:ngin
                ggrp = gin(gg);
                gvec = gvec + XS.ScatteringXS(cmat,grp,ggrp,k)*flux{ggrp,m}(cn);
            end
            src{g,m}(cn) = src{g,m}(cn) + M*gvec;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%