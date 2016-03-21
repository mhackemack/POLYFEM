%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve CFEM Transport Problem
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_cfem_transport(ndat, solvdat, mesh, DoF, FE, x, A)
global glob
% Setup Solution Space
% ------------------------------------------------------------------------------
flux_out = clear_flux_moments(x, DoF.TotalDoFs);
if ndat.Transport.HasReflectingBoundary || ndat.Transport.HasPeriodicBoundary
    ndat.Transport.ReflectingFluxesOld = ndat.Transport.ReflectingFluxes;
    ndat.Transport.PeriodicFluxesOld = ndat.Transport.PeriodicFluxes;
end
ndat.Transport.OutgoingCurrentsOld = ndat.Transport.OutgoingCurrents;
ndat.Transport.IncomingCurrentsOld = ndat.Transport.IncomingCurrents;
ndat = zero_partial_currents(ndat, mesh, DoF);
% Loop through Quadrature Directions
% ------------------------------------------------------------------------------
rev_str = [];
for m=1:ndat.Transport.NumberAngularDirections
    if glob.print_info
        msg = sprintf('   Calculating Flux for Angular Direction: %d',m);
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    
    [L, rhs] = get_global_matrices(m, ndat, mesh, DoF, FE, x);
    y = L\rhs;
    flux_out = add_to_flux(y, m, ndat, flux_out, DoF.TotalDoFs);
    if ndat.Transport.HasReflectingBoundary
        ndat = add_reflecting_angular_fluxes(ndat, mesh, DoF, m, y);
    end
    if ndat.Transport.HasPeriodicBoundary
        ndat = add_periodic_angular_fluxes(ndat, mesh, DoF, m, y);
    end
    % Add to the partial boundary currents
    ndat = compute_partial_boundary_currents(ndat, mesh, DoF, m, y);
end
% Perform DSA Update
% ------------------------------------------------------------------------------
if ndat.Transport.performDSA
    if glob.print_info
        msg = '   Performing DSA Correction Calculation.';
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    
    dx = cell(ndat.numberEnergyGroups,1);
    for g=1:ndat.numberEnergyGroups
        dx{g} = flux_out{g,1} - x{g,1};
    end
    [dx, A] = perform_cfem_DSA(ndat, solvdat, mesh, DoF, FE, dx, A);
    for g=1:ndat.numberEnergyGroups
        flux_out{g,1} = flux_out{g,1} + dx{g};
    end
end
% Outputs
% ------------------------------------------------------------------------------
nout = nargout;
varargout{1} = ndat;
varargout{2} = flux_out;
if nout > 2, varargout{3} = A; end
if nout > 3, varargout{4} = 0; end
if nout > 4, varargout{5} = 0; end
% Clear Command Line Text
if glob.print_info, fprintf(rev_str); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, rhs] = get_global_matrices(angNum, ndat, mesh, DoF, FE, x)
global glob
dim = mesh.Dimension;
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
angDir = ndat.Transport.AngularDirections(angNum,:);
m2d = ndat.Transport.moment_to_discrete(:,angNum);
% Allocate Memory
if ndof > glob.maxMatrix
    [L, rhs] = get_sparse_matrices(angNum, ndat, mesh, DoF, FE, x);
    return
else
    L = zeros(ndof,ndof);
    rhs = zeros(ndof,1);
end
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    F = FE.CellFunctionMatrix{tcell};
    GG = cell_dot(dim, angDir, FE.CellGradientMatrix{tcell});
    % Loop through energy groups
    for g=1:ng
        sg = cnodes + (g-1)*ndg;
        L(sg,sg) = L(sg,sg) + ndat.Transport.TotalXS(matID,g)*M + GG;
        % Add external source contribution
        if ndat.Transport.MMS
            gfunc = ndat.Transport.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx,angDir));
        else
            tvec = ndat.Transport.ExtSource(matID,g)*m2d(1) * F;
        end
        % Loop through cross energy groups
        for gg=1:ng
            fxs = ndat.Transport.FissSpec(matID,g)/ndat.keff*ndat.Transport.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg);
            tvec = tvec + (fxs*m2d(1))*M*x{gg,1}(cnodes);
            % Loop through scattering contributions from all flux moments
            for m=1:ndat.Transport.TotalFluxMoments
                k = ndat.Transport.MomentOrders(m,1) + 1;
                sxs = ndat.Transport.ScatteringXS(matID,gg,g,k)*m2d(m);
                tvec = tvec + sxs*M*x{gg,m}(cnodes);
            end
        end
        rhs(sg) = rhs(sg) + tvec;
    end
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    face = mesh.BoundaryFaces(f);
    fnorm = mesh.FaceNormal(face,:);
    if dot(angDir,fnorm) < 0
        flag = mesh.FaceID(face);
        tflag = ndat.Transport.BCFlags(flag);
        fnodes = DoF.getFaceCellNodes(face,1);
        for g=1:ng
            fcnodes = fnodes + (g-1)*ndg;
            L(fcnodes,:) = 0;
            for i=1:length(fcnodes)
                L(fcnodes(i),fcnodes(i)) = 1;
            end
            switch(tflag)
                case(glob.Vacuum)
                    rhs(fcnodes) = 0;
                case(glob.Reflecting)
                    opp_dir = ndat.Transport.ReflectingBoundaryAngles{face}(angNum);
                    opp_af = ndat.Transport.ReflectingFluxesOld{face}{opp_dir}(:,g);
                    rhs(fcnodes) = opp_af;
                case(glob.Periodic)
                    
                case(glob.IncidentIsotropic)
                    rhs(fcnodes) = ndat.Transport.BCVals{flag,g};
                case(glob.IncidentCurrent)
                    
                case(glob.IncidentBeam)
                    rhs(fcnodes) = ndat.Transport.BeamFluxes{f}(angNum, g);
                case(glob.Function)
                    fxn = DoF.NodeLocations(fnodes,:);
                    fvals = ndat.Transport.BCVals{flag,g}(fxn,angDir);
                    rhs(fcnodes) = fvals;
            end
        end
    end
end
% Set system matrix to sparse to make solve faster
if ~issparse(L), L = sparse(L); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, rhs] = get_sparse_matrices(angNum, ndat, mesh, DoF, FE, x)
global glob
dim = mesh.Dimension;
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
angDir = ndat.Transport.AngularDirections(angNum,:);
m2d = ndat.Transport.moment_to_discrete(:,angNum);
% Allocate Memory
rhs = zeros(ndof,1);
I = [];
J = [];
TMAT = [];
FDI = [];
FDJ = [];
FDMAT = [];
FRI = [];
FRJ = [];
FRMAT = [];

% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    ncnodes = length(cnodes);
    onesnodes = ones(ncnodes,1);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    F = FE.CellFunctionMatrix{tcell};
    GG = cell_dot(dim, angDir, FE.CellGradientMatrix{tcell});
    for g=1:ng
        sg = cnodes + (g-1)*ndg;
        rows = onesnodes*sg;
        cols = (onesnodes*sg)';
        tmat = ndat.Transport.TotalXS(matID,g)*M + GG;
        I = [I;rows(:)];
        J = [J;cols(:)];
        TMAT = [TMAT;tmat(:)];
        % Add external source contribution
        if ndat.Transport.MMS
            gfunc = ndat.Transport.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx,angDir));
        else
            tvec = ndat.Transport.ExtSource(matID,g)*m2d(1) * F;
        end
%         tvec = ndat.Transport.ExtSource(matID,g)*m2d(1) * F;
        % Loop through cross energy groups
        for gg=1:ng
            fxs = ndat.Transport.FissSpec(matID,g)/ndat.keff*ndat.Transport.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg);
            tvec = tvec + (fxs*m2d(1))*M*x{gg,1}(cnodes);
            % Loop through scattering contributions from all flux moments
            for m=1:ndat.Transport.TotalFluxMoments
                k = ndat.Transport.MomentOrders(m,1) + 1;
                sxs = ndat.Transport.ScatteringXS(matID,gg,g,k)*m2d(m);
                tvec = tvec + sxs*M*x{gg,m}(cnodes);
            end
        end
        rhs(sg) = rhs(sg) + tvec;
    end
end
L = sparse(I,J,TMAT,ndof,ndof);
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    face = mesh.BoundaryFaces(f);
    fnorm = mesh.FaceNormal(face,:);
    if dot(angDir,fnorm) < 0
        flag = mesh.FaceID(face);
        tflag = ndat.Transport.BCFlags(flag);
        fnodes = DoF.getFaceCellNodes(face,1);
        for g=1:ng
            fcnodes = fnodes + (g-1)*ndg;
            L(fcnodes,:) = 0;
            for i=1:length(fcnodes)
                L(fcnodes(i),fcnodes(i)) = 1;
            end
            switch(tflag)
                case(glob.Vacuum)
                    rhs(fcnodes) = 0;
                case(glob.Reflecting)
                    opp_dir = ndat.Transport.ReflectingBoundaryAngles{face}(angNum);
                    opp_af = ndat.Transport.ReflectingFluxesOld{face}{opp_dir}(:,g);
                    rhs(fcnodes) = opp_af;
                case(glob.Periodic)
                    
                case(glob.IncidentIsotropic)
                    rhs(fcnodes) = ndat.Transport.BCVals{flag,g};
                case(glob.IncidentCurrent)
                    
                case(glob.IncidentBeam)
                    rhs(fcnodes) = ndat.Transport.BeamFluxes{f}(angNum, g);
                case(glob.Function)
                    fxn = DoF.NodeLocations(fnodes,:);
                    fvals = ndat.Transport.BCVals{flag,g}(fxn,angDir);
                    rhs(fcnodes) = fvals;
            end
        end
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
function flux = add_to_flux(psi, m_dir, ndat, flux, ndofs)
for n=1:ndat.Transport.TotalFluxMoments
    d2m = ndat.Transport.discrete_to_moment(n,m_dir);
    for g=1:ndat.numberEnergyGroups
        flux{g,n} = flux{g,n} + psi((g-1)*ndofs+1:g*ndofs)*d2m;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = add_reflecting_angular_fluxes(ndat, mesh, DoF, m, x)
global glob
angDir = ndat.Transport.AngularDirections(m,:);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fflag = mesh.FaceID(ff);
    fdot = angDir*mesh.FaceNormal(ff,:)';
    if ndat.Transport.BCFlags(fflag) == glob.Reflecting && fdot > 0
        fnodes = DoF.FaceCellNodes{ff,1};
        % Loop through all indices and save values
        for g=1:ndat.numberEnergyGroups
            sg = fnodes + (g-1)*DoF.TotalDoFs;
            ndat.Transport.ReflectingFluxes{ff}{m}(:,g) = x(sg);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = add_periodic_angular_fluxes(ndat, mesh, DoF, m, x)
global glob
angDir = ndat.Transport.AngularDirections(m,:);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fflag = mesh.FaceID(ff);
    fdot = angDir*mesh.FaceNormal(ff,:)';
    if ndat.Transport.BCFlags(fflag) == glob.Periodic && fdot > 0
        pn = DoF.PeriodicFaceDoFs{ff};
        of = mesh.PeriodicOppositeFaces(ff);
        fnodes = DoF.FaceCellNodes{ff,1};
        for g=1:ndat.numberEnergyGroups
            sg = fnodes + (g-1)*DoF.TotalDoFs;
            ndat.Transport.PeriodicFluxes{of}{m}(:,g) = x(sg);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = compute_partial_boundary_currents(ndat, mesh, DoF, m, x)
global glob
angDir = ndat.Transport.AngularDirections(m,:);
wt = ndat.Transport.AngularWeights(m);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fdot = angDir*mesh.FaceNormal(ff,:)';
    afdot = abs(fdot);
    if abs(fdot) < glob.small, continue; end
    fnodes = DoF.FaceCellNodes{ff,1};
    for g=1:ndat.numberEnergyGroups
        sg = fnodes + (g-1)*DoF.TotalDoFs;
        if fdot > 0
            tv = ndat.Transport.OutgoingCurrents{ff};
            tv(:,g) = tv(:,g) + wt * afdot * x(sg);
            ndat.Transport.OutgoingCurrents{ff}(:,g) = tv(:,g);
        elseif fdot < 0
            tv = ndat.Transport.IncomingCurrents{ff};
            tv(:,g) = tv(:,g) + wt * afdot * x(sg);
            ndat.Transport.IncomingCurrents{ff}(:,g) = tv(:,g);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = zero_partial_currents(ndat, mesh, DoF)
ng = ndat.numberEnergyGroups;
ndat.Transport.OutgoingCurrents = cell(mesh.TotalFaces, 1);
ndat.Transport.IncomingCurrents = cell(mesh.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fdn = length(DoF.FaceCellNodes{ff,1});
    ndat.Transport.OutgoingCurrents{ff} = zeros(fdn,ng);
    ndat.Transport.IncomingCurrents{ff} = zeros(fdn,ng);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%