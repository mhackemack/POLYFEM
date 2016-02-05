%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          SI-upwind + MIP matrix functor
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
function P = func_build_SI_upwind_MIP(lam, input)
global glob
% Build Phase Matrix
if dim == size(lam,2); lam=lam'; end
node_locs = input.dof.NodeLocations;
PV = exp(1i*node_locs*lam);
PM = diag(PV);
% Retrieve terms
ndofs = input.dof.TotalDoFs;
ng = input.data.numberEnergyGroups;
g_offset = (1:ng)*ndofs - ndofs;
m2d = input.Quadrature.moment_to_discrete;
d2m = input.Quadrature.discrete_to_moment;
% Retrieve Matrices
T = func_mat_SI_upwind(lam, input);
A = func_mat_MIP(lam, input);
S0 = input.ScatteringMatrix;
I = eye(ndofs*ng);
% Build full matrix based on acceleration type
% P = T + (A\B)*(T - I);
if input.data.AccelType == glob.Accel_WGS_DSA
    S = zeros(ng*ndofs);
    for g=1:ng
        gdofs = (1:ndofs) + g_offset(g);
        for gg=1:ng
            ggdofs = (1:ndofs) + g_offset(gg);
            S(gdofs,ggdofs) = S0(:,:,g,gg)*PM;
        end
    end
elseif input.data.AccelType == glob.Accel_AGS_TG
    
elseif input.data.AccelType == glob.Accel_AGS_MTG
    
end
