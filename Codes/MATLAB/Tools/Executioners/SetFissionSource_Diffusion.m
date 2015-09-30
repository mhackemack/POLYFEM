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
function src = SetFissionSource_Diffusion(XS,groups,gin,flux,mesh,DoF,FE,keff)
% Get some preliminary information
% ------------------------------------------------------------------------------
if nargin < 9 || isempty(keff), keff = 1.0; end
ngroups = length(groups); ngin = length(gin);
fxs = XS.NuBar*XS.FissionXS/keff;
% Allocate memory space
% ------------------------------------------------------------------------------
src = cell(ngroups,1);
for g=1:ngroups
    src{g} = zeros(DoF.TotalDoFs,1);
end
% Exit with vectors of zeros if there is no inscattering
if ngin == 0 || isempty(fxs), return; end
% Loop through spatial cells and build source
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    cn   = DoF.ConnectivityArray{c}; ncn = length(cn);
    M    = FE.CellMassMatrix{c};
    % Double loop through energy groups
    for g=1:ngroups
        chi = XS.FissSpec(cmat,groups(g));
        gvec = zeros(ncn,1);
        for gg=1:ngin
            gvec = gvec + chi*fxs(cmat,ngin(gg))*flux{ngin(gg)}(cn);
        end
        src{g}(cn) = src{g}(cn) + M*gvec;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%