%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve DFEM Transport Problem
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_dfem_transport(ndat, solvdat, mesh, DoF, FE, x, A)
global glob
% Setup Solution Space
% --------------------
flux_out = clear_flux_moments(x, DoF.TotalDoFs);
if ndat.Transport.HasReflectingBoundary || ndat.Transport.HasPeriodicBoundary
    ndat.Transport.ReflectingFluxesOld = ndat.Transport.ReflectingFluxes;
    ndat.Transport.PeriodicFluxesOld = ndat.Transport.PeriodicFluxes;
end
ndat.Transport.OutgoingCurrentsOld = ndat.Transport.OutgoingCurrents;
ndat.Transport.IncomingCurrentsOld = ndat.Transport.IncomingCurrents;
ndat = zero_partial_currents(ndat, mesh, DoF);
% Loop through Quadrature Directions
% ----------------------------------
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
    % Switch based on DSA type
    if strcmp(ndat.Transport.DSAType, 'MIP')
        [dx, A] = perform_MIP_DSA(ndat, solvdat, mesh, DoF, FE, dx, A);
    elseif strcmp(ndat.Transport.DSAType, 'IP')
        [dx, A] = perform_IP_DSA(ndat, solvdat, mesh, DoF, FE, dx, A);
    end
    for g=1:ndat.numberEnergyGroups
        flux_out{g,1} = flux_out{g,1} + dx{g};
    end
end

% Outputs
nout = nargout;
varargout{1} = ndat;
varargout{2} = flux_out;
if nout == 3
    if exist('A', 'var')
        varargout{3} = A; 
    else
        varargout{3} = []; 
    end
end
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
tangDir = angDir';
angNorm = ndat.Transport.AngQuadNorm;
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
        tvec = ndat.Transport.ExtSource(matID,g)*m2d(1) * F;
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
% Loop through all faces
for f=1:mesh.TotalFaces
    fflag = mesh.FaceID(f);
    fdot = mesh.FaceNormal(f,:)*tangDir;
    if abs(fdot) < glob.small, continue; end
    % Perform Interior Face operations
    if fflag == 0
        % Fixup face normal and cell ordering - this is upwinding
        if fdot < 0
            fdot = -1*fdot;
            fcnodes1 = DoF.FaceCellNodes{f,1};
            fcnodes2 = DoF.FaceCellNodes{f,2};
            M = FE.FaceMassMatrix{f,1};
            MM = FE.FaceConformingMassMatrix{f,1};
        else
            fcnodes1 = DoF.FaceCellNodes{f,2};
            fcnodes2 = DoF.FaceCellNodes{f,1};
            M = FE.FaceMassMatrix{f,2};
            MM = FE.FaceConformingMassMatrix{f,2};
        end
        for g=1:ng
            gfnodes1 = fcnodes1 + (g-1)*ndg;
            gfnodes2 = fcnodes2 + (g-1)*ndg;
            L(gfnodes1,gfnodes1) = L(gfnodes1,gfnodes1) + fdot*M;
            L(gfnodes1,gfnodes2) = L(gfnodes1,gfnodes2) - fdot*MM;
        end
    else
        % Only add boundary contribution to incoming fluxes
        if fdot < 0
            tflag = ndat.Transport.BCFlags(fflag);
            M = FE.FaceMassMatrix{f,1};
            F = FE.FaceFunctionMatrix{f,1};
            fnodes = DoF.FaceCellNodes{f,1};
            for g=1:ng
                fcnodes = fnodes + (g-1)*ndg;
                L(fcnodes,fcnodes) = L(fcnodes,fcnodes) - M*fdot;
                switch(tflag)
                    case(glob.Vacuum)
                        % do nothing
                    case(glob.Reflecting)
                        opp_dir = ndat.Transport.ReflectingBoundaryAngles{f}(angNum);
                        opp_af = ndat.Transport.ReflectingFluxesOld{f}{opp_dir}(:,g);
                        rhs(fcnodes) = rhs(fcnodes) - fdot * M * opp_af;
                    case(glob.Periodic)
                        per_af = ndat.Transport.PeriodicFluxesOld{f}{angNum}(:,g);
                        rhs(fcnodes) = rhs(fcnodes) - fdot * M * per_af;
                    case(glob.IncidentIsotropic)
                        rhs(fcnodes) = rhs(fcnodes) - (fdot*ndat.Transport.BCVals{fflag,g}/angNorm)*F;
                    case(glob.IncidentCurrent)
                        
                    case(glob.IncidentBeam)
                        beam_val = ndat.Transport.BeamFluxes{f}(angNum, g);
                        rhs(fcnodes) = rhs(fcnodes) - (fdot*beam_val/angNorm)*F;
                end
            end
        end
    end
end
if ~issparse(L), L = sparse(L); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, rhs] = get_sparse_matrices(angNum, ndat, mesh, DoF, FE, x)
global glob
dim = mesh.Dimension;
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
angDir = ndat.Transport.AngularDirections(angNum,:);
tangDir = angDir';
angNorm = ndat.Transport.AngQuadNorm;
m2d = ndat.Transport.moment_to_discrete(:,angNum);
% Allocate Memory
% ---------------
rhs = zeros(ndof,1);
I = [];
J = [];
TMAT = [];

% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    onesnodes = ones(ncnodes,1);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    F = FE.CellFunctionMatrix{tcell};
    GG = cell_dot(dim, angDir, FE.CellGradientMatrix{tcell});
    for g=1:ng
        % Get Indexing Information
        % ------------------------
        sg = cnodes + (g-1)*ndg;
        rows = (onesnodes*sg)';
        cols = onesnodes*sg;
        % Apply Sparse Indexing
        % ---------------------
        tmat = ndat.Transport.TotalXS(matID,g)*M + GG;
        I = [I;rows(:)];
        J = [J;cols(:)];
        TMAT = [TMAT;tmat(:)];
        tvec = ndat.Transport.ExtSource(matID,g)*m2d(1) * F;
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
% Loop through all interior faces
% -------------------------------
for ff=1:mesh.TotalInteriorFaces
    f = mesh.InteriorFaces(ff);
    fflag = mesh.FaceID(f);
    % Test boundary type of face dot product
    fdot = angDir*mesh.FaceNormal(f,:)';
    if fflag ~= 0 || abs(fdot) < glob.small, continue, end
    % Fixup face normal and cell ordering
    if fdot < 0
        fdot = -1*fdot;
        fcnodes1 = DoF.FaceCellNodes{f,1};
        fcnodes2 = DoF.FaceCellNodes{f,2};
        M = FE.FaceMassMatrix{f,1};
        MM = FE.FaceConformingMassMatrix{f,1};
    else
        fcnodes1 = DoF.FaceCellNodes{f,2};
        fcnodes2 = DoF.FaceCellNodes{f,1};
        M = FE.FaceMassMatrix{f,2};
        MM = FE.FaceConformingMassMatrix{f,2};
    end
    onesnodes = ones(length(fcnodes1),1);
    for g=1:ng
        gfnodes1 = fcnodes1 + (g-1)*ndg;
        gfnodes2 = fcnodes2 + (g-1)*ndg;
        % Apply within cell term
        % ----------------------
        cols1 = onesnodes*gfnodes1; rows1 = (onesnodes*gfnodes1)'; tmat = fdot*M;
        I = [I;rows1(:)]; J = [J;cols1(:)]; TMAT = [TMAT;tmat(:)];
        % Apply cross-cell term
        % ---------------------
        cols2 = onesnodes*gfnodes2; tmat = -1.0*fdot*MM;
%         rows2 = (onesnodes*gfnodes2)';
        I = [I;rows1(:)]; J = [J;cols2(:)]; TMAT = [TMAT;tmat(:)];
    end
end
% Create sparse system matrix
% ---------------------------
L = sparse(I,J,TMAT,ndof,ndof);
%
% Loop through all boundary faces
% -------------------------------
for ff=1:mesh.TotalBoundaryFaces
    f = mesh.BoundaryFaces(ff);
    fflag = mesh.FaceID(f);
    fdot = mesh.FaceNormal(f,:)*tangDir;
    % Only add boundary contribution to incoming fluxes
    if fdot < 0 && abs(fdot) > glob.small
        tflag = ndat.Transport.BCFlags(fflag);
        M = FE.FaceMassMatrix{f,1};
        F = FE.FaceFunctionMatrix{f,1};
        fnodes = DoF.FaceCellNodes{f,1};
        for g=1:ng
            fcnodes = fnodes + (g-1)*ndg;
            L(fcnodes,fcnodes) = L(fcnodes,fcnodes) - M*fdot;
            switch(tflag)
                case(glob.Vacuum)
                    % do nothing
                case(glob.Reflecting)
                    opp_dir = ndat.Transport.ReflectingBoundaryAngles{f}(angNum);
                    opp_af = ndat.Transport.ReflectingFluxesOld{f}{opp_dir}(:,g);
                    rhs(fcnodes) = rhs(fcnodes) - fdot * M * opp_af;
                case(glob.Periodic)
                    per_af = ndat.Transport.PeriodicFluxesOld{f}{angNum}(:,g);
                    rhs(fcnodes) = rhs(fcnodes) - fdot * M * per_af;
                case(glob.IncidentIsotropic)
                    rhs(fcnodes) = rhs(fcnodes) - (fdot*ndat.Transport.BCVals{fflag,g}/angNorm)*F;
                case(glob.IncidentCurrent)
                    
                case(glob.IncidentBeam)
                    beam_val = ndat.Transport.BeamFluxes{f}(angNum, g);
                    rhs(fcnodes) = rhs(fcnodes) - (fdot*beam_val/angNorm)*F;
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