%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execute Sweep Chunk
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    Invert Transport Operator by Sweeping - upwind scheme is
%                   strongly enforced.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flux = exec_func_transport_sweep(ndat, mesh, DoF, FE, x, m, groups)
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
CellOrder = sweep.CellSweepOrder{m}; ncells = length(CellOrder);
USFaces = sweep.UpstreamFaces{m};
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
        GG = cell_dot(dim, adir, G);
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
            A(:,:,q,g) = txs(cmat,grp)*M + GG;
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
        nfnodes = length(DoF.FaceCellNodes{f,1});
        fid = mesh.FaceID(f);
        fcells = mesh.FaceCells(f,:);
        fnorm = mesh.FaceNormal(f,:)';
        % Interior Faces
        if fid == 0
            if fcells(1) == c
                fn1 = DoF.FaceCellNodeNumbering{f,1};
                fn2 = (DoF.FaceCellNodes{f,2});
                fn2 = fn2(nfnodes:-1:1);
                M = FE.FaceMassMatrix{f,1};
            else
                fn1 = DoF.FaceCellNodeNumbering{f,2};
                fn2 = (DoF.FaceCellNodes{f,1});
                fn2 = fn2(nfnodes:-1:1);
                M = FE.FaceMassMatrix{f,2};
                fnorm = -fnorm;
            end
            % Loop through angles
            for q=1:na
                tq = angs(q);
                adir = angdirs(tq,:);
%                 fdot = adir*fnorm;
                MMM = adir*fnorm*M;
                % Loop through energy groups
                for g=1:ng
%                     psi = flux(fn2,q,g);
                    A(fn1,fn1,q,g) = A(fn1,fn1,q,g) - MMM;
                    b(fn1,q,g) = b(fn1,q,g) - MMM*flux(fn2,q,g);
                end
            end
        % Boundary Faces
        else
            tflag = ndat.Transport.BCFlags(fid);
            M = FE.FaceMassMatrix{f,1};
            fnodes = DoF.FaceCellNodes{f,1};
            fn = DoF.FaceCellNodeNumbering{f,1}; 
            ofn = ones(length(fn), 1); zfn = zeros(length(fn), 1);
            % Loop through angles
            for q=1:na
                tq = angs(q);
                adir = angdirs(tq,:);
                MM = adir*fnorm*M;
                % Loop through energy groups
                for g=1:ng
                    grp = groups(g);
                    A(fn,fn,q,g) = A(fn,fn,q,g) - MM;
                    switch(tflag)
                        case(glob.Vacuum)
                            psi = zfn;
                        case(glob.Reflecting)
                            opp_dir = ndat.Transport.ReflectingBoundaryAngles{f}(tq);
                            psi = ndat.Transport.ReflectingFluxes{f}{opp_dir}(:,g);
                        case(glob.IncidentIsotropic)
                            val = ndat.Transport.BCVals{fid,grp}/angNorm;
                            psi = val*ofn;
                        case(glob.IncidentBeam)
                            beam_val = ndat.Transport.BeamFluxes{f}(tq, g);
                            psi = beam_val*ofn;
                        case(glob.Function)
                            fxn = DoF.NodeLocations(fnodes,:);
                            psi = ndat.Transport.BCVals{fid,grp}(fxn,adir);
                    end
                    b(fn,q,g) = b(fn,q,g) - MM*psi;
                end
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
function out = cell_dot(dim, vec1, vec2)
if dim == 1
    out = vec1*vec2{1};
elseif dim == 2
    out = vec1(1)*vec2{1} + vec1(2)*vec2{2};
else
    out = vec1(1)*vec2{1} + vec1(2)*vec2{2} + vec1(3)*vec2{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
