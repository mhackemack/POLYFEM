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
function E = func_build_SI_upwind_MIP(lam, input)
global glob
% Build Phase Matrix
if input.mesh.Dimension == size(lam,2); lam=lam'; end
node_locs = input.dof.NodeLocations;
PV = exp(1i*node_locs*lam);
PM = diag(PV);
% Retrieve terms
ng = input.data.numberEnergyGroups;
ndofs = input.dof.TotalDoFs;
zdofs = (1:ndofs);
g_offset = (1:ng)*ndofs - ndofs;
m2d = input.Quadrature.moment_to_discrete;
d2m = input.Quadrature.discrete_to_moment;
% Retrieve Matrices
L = func_mat_SI_upwind(lam, input);
A = func_mat_MIP(lam, input);
S0 = input.ScatteringMatrix;
P = input.ProjectionMatrix;
I = speye(ndofs*ng);
% Build full matrix based on acceleration type
if input.data.AccelType == glob.Accel_WGS_DSA
    % Build Scattering matrix
    S = zeros(ng*ndofs);
    for g=1:ng
        gdofs = zdofs + g_offset(g);
        for gg=1:ng
            ggdofs = zdofs + g_offset(gg);
            S(gdofs,ggdofs) = S0(:,:,g,gg)*PM;
        end
    end
    % Build transport operators
    T = zeros(ndofs*ng);
    for q=1:input.Quadrature.NumberAngularDirections
        T = T + d2m(1,q)*( L{q}\( m2d(1,q)*S ));
    end
    TT = T - I;
    % Build final matrix based on energy group numbers
    if ng == 1
        E = T + (A\S)*TT;
    else
        MDSA = zeros(ndofs,ng*ndofs); SS = S*TT;
        for g=1:ng
            gdofs = zdofs + g_offset(g);
            for gg=(g+1):ng
            end
            MDSA = MDSA + SS(gdofs,:);
        end
        E = T + P*(A\MDSA);
    end
elseif input.data.AccelType == glob.Accel_WGS_INNER_DSA
    Sl = zeros(ng*ndofs);
    Sd = zeros(ng*ndofs);
    Su = zeros(ng*ndofs);
    for g=1:ng
        gdofs = zdofs + g_offset(g);
        for gg=1:ng
            ggdofs = zdofs + g_offset(gg);
            tS0 = S0(:,:,g,gg)*PM;
            if gg < g
                Sl(gdofs,ggdofs) = tS0;
            elseif gg==g
                Sd(gdofs,ggdofs) = tS0;
            elseif gg > g
                Su(gdofs,ggdofs) = tS0;
            end
        end
    end
    Sl = sparse(Sl); Sd = sparse(Sd); Su = sparse(Su); S = Sl + Sd + Su;
    % Build tranpsort matrices
    T = zeros(ndofs*ng);
    for q=1:input.Quadrature.NumberAngularDirections
        T = T + d2m(1,q)*( L{q}\( m2d(1,q)*S ));
    end
    TT = T - I;
    MDSA = zeros(ndofs,ng*ndofs); SS = (Sl+Su)*TT;
    
elseif input.data.AccelType == glob.Accel_AGS_TG
    % Build scattering matrices
    Sd = zeros(ng*ndofs);
    Su = zeros(ng*ndofs);
    for g=1:ng
        gdofs = zdofs + g_offset(g);
        for gg=1:ng
            ggdofs = zdofs + g_offset(gg);
            tS0 = S0(:,:,g,gg)*PM;
            if gg <= g
                Sd(gdofs,ggdofs) = tS0;
            elseif gg > g
                Su(gdofs,ggdofs) = tS0;
            end
        end
    end
    Sd = sparse(Sd); Su = sparse(Su);
    % Build tranpsort matrices
    B = zeros(ng*ndofs); C = zeros(ng*ndofs);
    for q=1:input.Quadrature.NumberAngularDirections
        tL = sparse(L{q}\I);
        B = B + d2m(1,q)*( tL*( m2d(1,q)*Sd ));
        C = C + d2m(1,q)*( tL*( m2d(1,q)*Su ));
    end
    TC = (I-B)\C; TCI = TC - I;
    % Build restriction/projection operators
    MDSA = zeros(ndofs,ng*ndofs); SS = Su*TCI;
    for g=1:ng
        gdofs = zdofs + g_offset(g);
        MDSA = MDSA + SS(gdofs,:);
    end
    E = TC + P*(A\MDSA);
elseif input.data.AccelType == glob.Accel_AGS_MTG
    % Build scattering matrices
    Sd = zeros(ng*ndofs);
    Su = zeros(ng*ndofs);
    for g=1:ng
        gdofs = zdofs + g_offset(g);
        for gg=1:ng
            ggdofs = zdofs + g_offset(gg);
            tS0 = S0(:,:,g,gg)*PM;
            if gg < g
                Sd(gdofs,ggdofs) = tS0;
            elseif gg >= g
                Su(gdofs,ggdofs) = tS0;
            end
        end
    end
    Sd = sparse(Sd); Su = sparse(Su);
    % Build tranpsort matrices
    B = zeros(ng*ndofs); C = zeros(ng*ndofs);
    for q=1:input.Quadrature.NumberAngularDirections
        tL = sparse(L{q}\I);
        B = B + d2m(1,q)*( tL*( m2d(1,q)*Sd ));
        C = C + d2m(1,q)*( tL*( m2d(1,q)*Su ));
    end
    TC = (I-B)\C; TCI = TC - I;
    % Build restriction/projection operators
    MDSA = zeros(ndofs,ng*ndofs); SS = Su*TCI;
    for g=1:ng
        gdofs = zdofs + g_offset(g);
        MDSA = MDSA + SS(gdofs,:);
    end
    E = TC + P*(A\MDSA);
end
