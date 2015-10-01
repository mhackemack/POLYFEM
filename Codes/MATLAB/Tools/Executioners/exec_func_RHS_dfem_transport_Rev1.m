%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execution Functor - RHS DFEM Transport
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = exec_func_RHS_dfem_transport_Rev1(data,qid,angs,groups,mesh,DoF,FE)
% Process Input Space
% ------------------------------------------------------------------------------
ndof = DoF.TotalDoFs;
ng = length(groups);
na = length(angs);
ntot = ndof * ng * na;
mquad = data.Quadrature(qid);
angdirs = mquad.AngularDirections';
m2d = mquad.moment_to_discrete;
% Build rhs vector and some stride information
% ------------------------------------------------------------------------------
rhs = zeros(ntot, 1);
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
% Loop through angles and groups to build volumetric right hand side terms
% ------------------------------------------------------------------------------
for q=1:na
    tq = angs(q);
    for g=1:ng
        grp = groups(g);
        qg_off = q_offset(q) + g_offset(g);
        % Apply Fission Terms - 0th order only
        rhs(1:ndof+qg_off) = rhs(1:ndof+qg_off) + data.Sources.FissionSource{grp,1}*m2d(1,tq);
        % Apply External Source Terms - 0th order only
        rhs(1:ndof+qg_off) = rhs(1:ndof+qg_off) + data.Sources.ExtSource{grp,1}*m2d(1,tq);
        for m=1:mquad.TotalFluxMoments
            % Apply InScattering Terms
            rhs(1:ndof+qg_off) = rhs(1:ndof+qg_off) + data.Sources.InScatteringSource{grp,m}*m2d(m,tq);
            % Apply WithinScattering Terms
            rhs(1:ndof+qg_off) = rhs(1:ndof+qg_off) + data.Sources.WithinScatteringSource{grp,m}*m2d(m,tq);
        end
    end
end 
% Loop through boundary faces and set incident flux contributions
% ------------------------------------------------------------------------------
for ff=1:mesh.TotalBoundaryFaces
    f = mesh.BoundaryFaces(ff);
    inc_flux = data.Fluxes.IncomingBoundaryFlux{f};
    fnorm = mesh.FaceNormal(f,:);
    fnodes = DoF.FaceCellNodes{f,1};
    M = FE.FaceMassMatrix{f,1};
    % Loop through angles
    for q=1:na
        tq = angs(q);
        fdot = fnorm*angdirs(:,tq);
        if fdot > 0, continue; end
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
            cnqg = fnodes + g_offset(g) + q_offset(q);
            rhs(cnqg) = rhs(cnqg) - fdot*M*inc_flux{tq,grp};
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%