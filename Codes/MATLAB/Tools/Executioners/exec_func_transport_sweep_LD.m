%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execute LD Sweep Chunk
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    Invert Transport Operator by Sweeping - upwind scheme is
%                   strongly enforced.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flux = exec_func_transport_sweep_LD(ndat, mesh, DoF, FE, x, m, groups)
% Process Input Space
% ------------------------------------------------------------------------------
global glob
dim = mesh.Dimension;
ndof = DoF.TotalDoFs;
ng = length(groups);
angs = ndat.Transport.AngleSets{m}; na = length(angs);
angdirs = ndat.Transport.AngularDirections;
angNorm = ndat.Transport.AngQuadNorm;
m2d = ndat.Transport.moment_to_discrete;
Kn = ndat.Transport.MomentOrders(:,1) + 1;
% Get Sweep Information
% ------------------------------------------------------------------------------
sweep = ndat.Transport.Sweeping;
CellOrder = sweep.CellSweepOrder{m};
USFaces = sweep.UpstreamFaces{m};
USCells = sweep.UpstreamCells{m};
DSFaces = sweep.DownstreamFaces{m};
DSCells = sweep.DownstreamCells{m};
% Allocate Memory Space
% ------------------------------------------------------------------------------
flux = zeros(ndof, na, ng);
txs = ndat.Transport.TotalXS;
sxs = ndat.Transport.ScatteringXS;
Qext = ndat.Transport.ExtSource;
keff = ndat.keff;
fxs = ndat.Transport.FissionXS;
chi = ndat.Transport.FissSpec;
nu = ndat.Transport.NuBar;
% Update Sweep Plot if Boolean Set
% ------------------------------------------------------------------------------
if ndat.Transport.visualizeSweeping
    initialize_sweep_plot(mesh);
end
% Loop through Cells in Sweep Chunk
% ------------------------------------------------------------------------------
for cc=1:mesh.TotalCells
    c = CellOrder(cc);
    cmat = mesh.MatID(c);
    cnodes = DoF.ConnectivityArray{c}; tcnodes = cnodes'; ncnodes = length(cnodes);
    M = FE.CellMassMatrix{c};
    F = FE.CellIntegralVector{c};
    G = FE.CellGradientMatrix{c};
    A = zeros(ncnodes,ncnodes,na,ng); b = zeros(ncnodes,na,ng);
    % Loop through angles
    % --------------------------------------------------------------------------
    for q=1:na
        tq = angs(q);
        adir = angdirs(tq,:);
%         GG = cell_dot(dim, adir, G);
        GG = cell_dot(dim, adir, G)';
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
%             A(:,:,q,g) = txs(cmat,grp)*M + GG;
            A(:,:,q,g) = txs(cmat,grp)*M - GG;
            if ndat.Transport.MMS
                gfunc = ndat.Transport.ExtSource{g};
                qx = FE.CellQuadNodes{c};
                qw = FE.CellQuadWeights{c};
                cb = FE.CellBasisValues{c};
                tvec = cb'*(qw.*gfunc(qx,adir));
            else
                tvec = Qext(cmat,grp)*m2d(1,tq) * F;
            end
            % Loop through energy groups again
            for gg=1:ng
                ggrp = groups(gg);
                % apply fission term
                tfxs = chi(cmat,grp)/keff*fxs(cmat,ggrp)*nu(cmat,ggrp);
                tvec = tvec + tfxs*M*x{ggrp,1}(cnodes);
                % Add scattering contribution
                for m=1:ndat.Transport.TotalFluxMoments
                    k = Kn(m);
                    tsxs = sxs(cmat,ggrp,grp,k)*m2d(m,tq);
                    tvec = tvec + tsxs*M*x{ggrp,m}(cnodes);
                end
            end
            b(:,q,g) = tvec;
        end
    end
    % Loop through upstream faces
    % --------------------------------------------------------------------------
    cfaces = USFaces{c};
    for ff=1:length(cfaces);
        f = cfaces(ff);
        fid = mesh.FaceID(f);
        fcells = mesh.FaceCells(f,:);
        fnorm = mesh.FaceNormal(f,:)';
        % Interior Faces
        if fid == 0
            if fcells(1) == c
%                 fn2 = DoF.FaceCellNodes{f,2};
                invec = FE.FaceIntegralVector{f,1};
                upvec = FE.FaceIntegralVector{f,2};
                updof = DoF.FaceCellNodes{f,2};
            else
%                 fn2 = DoF.FaceCellNodes{f,1};
                invec = FE.FaceIntegralVector{f,2};
                upvec = FE.FaceIntegralVector{f,1};
                updof = DoF.FaceCellNodes{f,1};
                fnorm = -1*fnorm;
            end
            M = invec*upvec';
            % Loop through angles
            for q=1:na
                tq = angs(q);
                adir = angdirs(tq,:);
                MMM = adir*fnorm*M;
                % Loop through energy groups
                for g=1:ng
                    upflux = flux(updof,q,g);
                    b(:,q,g) = b(:,q,g) - MMM*upflux;
                end
            end
        % Boundary Faces
        else
            tflag = ndat.Transport.BCFlags(fid);
            M = FE.FaceMassMatrix{f,1};
            F = FE.FaceIntegralVector{f,1};
            % Loop through angles
            for q=1:na
                tq = angs(q);
                adir = angdirs(tq,:);
                fdot = adir*fnorm;
                for g=1:ng
                    grp = groups(g);
                    switch(tflag)
                        case(glob.Vacuum)
                            % do nothing
                        case(glob.Reflecting)
                            opp_dir = ndat.Transport.ReflectingBoundaryAngles{f}(tq);
                            psi = ndat.Transport.ReflectingFluxes{f}{opp_dir}(:,g);
                            MM = fdot*M;
                            b(cnodes,q,g) = b(cnodes,q,g) - MM*psi;
                        case(glob.IncidentIsotropic)
                            val = ndat.Transport.BCVals{fid,grp}/angNorm;
                            b(cnodes,q,g) = b(cnodes,q,g) - fdot*F*val;
                    end
                end
            end
        end
    end
    % Loop through downstream faces
    % --------------------------------------------------------------------------
    cfaces = DSFaces{c};
    for ff=1:length(cfaces);
        f = cfaces(ff);
        fcells = mesh.FaceCells(f,:);
        fnorm = mesh.FaceNormal(f,:)';
        if fcells(1) == c
            M = FE.FaceMassMatrix{f,1};
        elseif fcells(2) == c
            fnorm = -1*fnorm;
            M = FE.FaceMassMatrix{f,2};
        end
        % Loop through angles
        for q=1:na
            tq = angs(q);
            adir = angdirs(tq,:);
            MMM = adir*fnorm*M;
            % Loop through energy groups
            for g=1:ng
                A(:,:,q,g) = A(:,:,q,g) + MMM;
            end
        end
    end
    % Loop through angles/groups and calculate new angular fluxes
    % --------------------------------------------------------------------------
    for q=1:na
        for g=1:ng
            flux(tcnodes,q,g) = A(:,:,q,g)\b(:,q,g);
        end
    end
    % Update Sweep Plot if Boolean Set
    % --------------------------------------------------------------------------
    if ndat.Transport.visualizeSweeping
        update_sweep_plot(mesh, c);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Listing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = cell_dot(dim, vec1, vec2)
if dim == 1
    out = vec1*vec2{1};
elseif dim == 2
    out = vec1(1)*vec2{1} + vec1(2)*vec2{2};
else
    out = vec1(1)*vec2{1} + vec1(2)*vec2{2} + vec1(3)*vec2{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%