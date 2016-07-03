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
for cc=1:ncells
    c = CellOrder(cc);
    cmat = mesh.MatID(c);
    cnodes = DoF.ConnectivityArray{c}; tcnodes = cnodes'; ncnodes = length(cnodes);
    M = FE.CellMassMatrix{c};
    F = FE.CellFunctionMatrix{c};
    G = FE.CellGradientMatrix{c};
    A = zeros(ncnodes,ncnodes,na,ng); b = zeros(ncnodes,na,ng);
    % Loop through angles
    % --------------------------------------------------------------------------
    for q=1:na
        tq = angs(q);
        adir = angdirs(tq,:);
        GG = cell_dot(dim, adir, G)';
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
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
    % Loop through faces by streaming type
    % --------------------------------------------------------------------------
    % Loop through upstream faces
    cfaces = USFaces{c};
    for ff=1:length(cfaces);
        f = cfaces(ff);
        fid = mesh.FaceID(f);
        fcells = mesh.FaceCells(f,:);
        fnorm = mesh.FaceNormal(f,:)';
        % Interior Faces
        if fid == 0
            if fcells(1) == c
                fn2 = DoF.FaceCellNodes{f,2};
            else
                fn2 = DoF.FaceCellNodes{f,1};
                fnorm = -1*fnorm;
            end
        % Boundary Faces
        else
            
        end
        % Loop through angles
        for q=1:na
            tq = angs(q);
            adir = angdirs(tq,:);
            
            % Loop through energy groups
            for g=1:ng
                
            end
        end
    end
    % Loop through downstream faces
    cfaces = DSFaces{c};
    for ff=1:length(cfaces);
        f = cfaces(ff);
        fnorm = mesh.FaceNormal(f,:)';
        M = FE.FaceMassMatrix{f,1};
        % Loop through angles
        for q=1:na
            tq = angs(q);
            adir = angdirs(tq,:);
            MMM = adir*fnorm*M;
            A(fn1,fn1,q,g) = A(fn1,fn1,q,g) + MMM;
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