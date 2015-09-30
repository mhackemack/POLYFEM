%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set Fission Diffusion Source
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function src = SetFissionSource_Diffusion(data,xsid,groups,flux,mesh,DoF,FE,keff)
% Get some preliminary information
% ------------------------------------------------------------------------------
if nargin < 8 || isempty(keff), keff = 1.0; end
ngroups = length(groups);
ntotg = data.Groups.NumberEnergyGroups;
XS = data.XS(xsid);
fxs = XS.NuBar*XS.FissionXS/keff;
% Allocate memory space
% ------------------------------------------------------------------------------
src = cell(ngroups,1);
% Loop through spatial cells and build source
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    cn   = DoF.ConnectivityArray{c}; ncn = length(cn);
    M    = FE.CellMassMatrix{c};
    % Double loop through energy groups
    for g=1:ngroups
        grp = groups(g);
        chi = XS.FissSpec(cmat,grp);
        gvec = zeros(ncn,1);
        for gg=1:ntotg
            gvec = gvec + chi*fxs(cmat,gg)*M*flux{gg}(cn);
        end
        src{g}(cn) = src{g}(cn) + gvec;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%