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
function rhs = exec_func_RHS_dfem_transport_Rev1(data)
% Process Input Space
% -------------------
ndof = DoF.TotalDoFs;
ng = length(groups);
na = length(angs);
ntot = ndof * ng * na;
angdirs = mquad.AngularDirections';
m2d = mquad.moment_to_discrete;
% Build rhs vector and some stride information
% --------------------------------------------
rhs = zeros(ntot, 1);
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
% Loop through cells
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    cnodes = DoF.ConnectivityArray{c};
    % Loop through angles
    for q=1:na
        tq = angs(q);
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
            cnqg = cnodes + g_offset(g) + q_offset(q);
            
        end
    end
end
% Loop through boundary faces
for ff=1:mesh.TotalBoundaryFaces
    f = mesh.BoundaryFaces(ff);
    inc_flux = data.Fluxes.IncomingBoundaryFlux{f};
    fnorm = mesh.FaceNormal(f,:);
    fnodes = DoF.FaceCellNodes{f,1};
    % Loop through angles
    for q=1:na
        tq = angs(q);
        fdot = fnorm*angdirs(:,tq);
        if fdot > 0, continue; end
        M = FE.FaceMassMatrix{f,1};
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
            cnqg = fnodes + g_offset(g) + q_offset(q);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%