%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execution Functor - RHS DFEM Transport (upwind)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = exec_func_RHS_dfem_transport_Rev1(x, data, mesh, DoF, FE, angs, groups)
% Process Input Space
% -------------------
global glob
ndof = DoF.TotalDoFs;
ng = length(groups);
na = length(angs);
ntot = ndof * ng * na;
keff = data.keff;
xs = data.XS;
mquad = data.Quadrature;
fluxes = data.Fluxes;
angdirs = mquad.AngularDirections';
angNorm = mquad.AngQuadNorm;
% Build rhs vector and some stride information
% --------------------------------------------
rhs = zeros(ntot, 1);
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
% Loop through cells
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    cnodes = DoF.ConnectivityArray{c};
    M = FE.CellMassMatrix{c};
    F = FE.CellFunctionMatrix{c};
    % Loop through angles
    for q=1:na
        tq = angs(q);
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
            cnqg = cnodes + g_offset(g) + q_offset(q);
            % Add external source contribution
            if data.MMS
                gfunc = xs.ExtSource{grp};
                qx = FE.CellQuadNodes{c};
                qw = FE.CellQuadWeights{c};
                cb = FE.CellBasisValues{c};
                tvec = cb'*(qw.*gfunc(qx,angdirs(:,tq)));
            else
                tvec = xs.ExtSource(cmat,grp)*mquad.moment_to_discrete(1,tq)*F;
            end
            % COMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
            % COMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
            % Loop through energy groups again
            for gg=1:ng
                ggrp = groups(gg);
                % apply fission term
                fxs = xs.FissSpec(cmat,grp)/keff*xs.FissionXS(cmat,ggrp)*xs.NuBar(cmat,ggrp);
                tvec = tvec + fxs*M*x{ggrp,1}(cnodes);
                % Add scattering contribution
                for m=1:fluxes.TotalFluxMoments
                    k = fluxes.MomentOrders(m,1) + 1;
                    sxs = xs.ScatteringXS(cmat,ggrp,grp,k)*mquad.moment_to_discrete(m,tq);
                    tvec = tvec + sxs*M*x{ggrp,m}(cnodes);
                end
            end
            % COMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
            % COMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
            % Apply local matrix contribution
            rhs(cnqg) = rhs(cnqg) + tvec;
        end
    end
end
% Loop through boundary faces
for ff=1:mesh.TotalBoundaryFaces
    f = mesh.BoundaryFaces(ff);
    fflag = mesh.FaceID(f);
    tflag = xs.BCFlags(fflag);
    fnorm = mesh.FaceNormal(f,:);
    fnodes = DoF.FaceCellNodes{f,1};
    % Loop through angles
    for q=1:na
        tq = angs(q);
        fdot = fnorm*angdirs(:,tq);
        if fdot > 0, continue; end
        M = FE.FaceMassMatrix{f,1};
        F = FE.FaceFunctionMatrix{f,1};
        % Loop through energy groups
        for g=1:ng
            grp = groups(g);
            cnqg = fnodes + g_offset(g) + q_offset(q);
            switch(tflag)
                case(glob.Reflecting)
                    opp_dir = data.Fluxes.ReflectingBoundaryAngles{f}(tq);
                    opp_af = data.Fluxes.ReflectingFluxes{f}{opp_dir}(:,grp);
                    rhs(cnqg) = rhs(cnqg) - fdot * M * opp_af;
                case(glob.Periodic)
                    per_af = data.Fluxes.PeriodicFluxesOld{f}{tq}(:,grp);
                    rhs(cnqg) = rhs(cnqg) - fdot * M * per_af;
                case(glob.IncidentIsotropic)
                    rhs(cnqg) = rhs(cnqg) - (fdot*xs.BCVals(fflag,grp)/angNorm)*F;
                case(glob.IncidentCurrent)

                case(glob.IncidentBeam)
                    beam_val = data.Fluxes.BeamFluxes{f}(tq, grp);
                    rhs(cnqg) = rhs(cnqg) - (fdot*beam_val/angNorm)*F;
                case(glob.Function)
                    fxn = DoF.NodeLocations(fnodes,:);
                    fvals = xs.BCVals{fflag,grp}(fxn,angdirs(:,tq));
                    rhs(cnqg) = rhs(cnqg) - fdot * M * fvals;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%