%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve DFEM Diffusion Problem
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_dfem_diffusion(ndat, solvdat, mesh, DoF, FE, x, A)
global glob
ndg = DoF.TotalDoFs;
ndof = ndat.numberEnergyGroups * ndg;
if nargin < 5 || isempty(x)
    x = ones(ndof,1);
else
    x = cell_to_vector(x, DoF);
end
if length(x) > glob.maxSparse
    rhs = get_rhs(x, ndat, mesh, DoF, FE);
    [x,flag,relres,iter] = pcg(@(x) get_Ax(x, ndat, mesh, DoF, FE),rhs,solvdat.relativeTolerance,solvdat.maxIterations,[],[],x);
else
    if nargin < 6 || ~exist('A','var')
        [A,rhs] = get_global_matrices(x, ndat, mesh, DoF, FE);
    else
        if isempty(A)
            [A,rhs] = get_global_matrices(x, ndat, mesh, DoF, FE);
        else
            rhs = get_rhs(x, ndat, mesh, DoF, FE);
        end
    end
    % Compute Diffusion Solution
    x = A\rhs;
end
x = vector_to_cell(x,DoF);
% Outputs
nout = nargout;
varargout{1} = ndat;
varargout{2} = x;
if nout == 3
    if exist('A', 'var')
        varargout{3} = A; 
    else
        varargout{3} = []; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,rhs] = get_global_matrices(x, ndat, mesh, DoF, FE)
global glob
dim = mesh.Dimension;
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
C_IP = ndat.IP_Constant;
% Allocate Memory
if ndof > glob.maxMatrix
    [L, rhs] = get_sparse_matrices(x, ndat, mesh, DoF, FE);
    return
else
    L = zeros(ndof,ndof);
    rhs = zeros(ndof,1);
end
rev_str = [];
% Loop through cells
for tcell=1:mesh.TotalCells
    % Print Current Cell Information
%     if glob.print_info
%         msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        L(sg,sg) = (L(sg,sg) + ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M);
        if ndat.Diffusion.MMS
            gfunc = ndat.Diffusion.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx));
        else
            tvec = ndat.Diffusion.ExtSource(matID,g)*M*ones(ncnodes,1);
        end
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (gg-1)*ndg;
            xsec = ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg) + ndat.Diffusion.ScatteringXS(matID,gg,g);
            tvec = tvec + xsec*M*x(sgg);
        end
    end
    rhs(sg) = rhs(sg) + tvec;
end
% Loop through faces
for f=1:mesh.TotalFaces
    % Print Current Face Information
%     if glob.print_info
%         msg = sprintf('   -> Building Face: %d of %d',f,mesh.TotalFaces);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    fflag = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:); tfnorm = fnorm';
    % Interior Face
    if fflag == 0
        matids = mesh.MatID(fcells);
        D = ndat.Diffusion.DiffXS(matids,:)';
        h = mesh.OrthogonalProjection(f,:);
        fcnodes1 = DoF.FaceCellNodes{f,1}; fcnodes2 = DoF.FaceCellNodes{f,2};
        cnodes1 = DoF.ConnectivityArray{fcells(1)}; cnodes2 = DoF.ConnectivityArray{fcells(2)};
        M1 = FE.FaceMassMatrix{f,1}; M2 = FE.FaceMassMatrix{f,2};
        MM1 = FE.FaceConformingMassMatrix{f,1}; MM2 = FE.FaceConformingMassMatrix{f,2};
        G1 = FE.FaceGradientMatrix{f,1}; G1 = cell_dot(dim,fnorm,G1);
        G2 = FE.FaceGradientMatrix{f,2}; G2 = cell_dot(dim,fnorm,G2);
        CG1 = FE.FaceCouplingGradientMatrix{f,1}; CG1 = cell_dot(dim,fnorm,CG1);
        CG2 = FE.FaceCouplingGradientMatrix{f,2}; CG2 = cell_dot(dim,fnorm,CG2);
        % Apply Interior Terms
        for g=1:ndat.numberEnergyGroups
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D(g,:), h, fflag);
            gfnodes1 = fcnodes1 + (g-1)*ndg;
            gfnodes2 = fcnodes2 + (g-1)*ndg;
            gcnodes1 =  cnodes1 + (g-1)*ndg;
            gcnodes2 =  cnodes2 + (g-1)*ndg;
            % Mass Terms
            % ------------------------------------------------------------------
            % (-,-)
            L(gfnodes1,gfnodes1) = L(gfnodes1,gfnodes1) + kp*M1;
            % (+,+)
            L(gfnodes2,gfnodes2) = L(gfnodes2,gfnodes2) + kp*M2;
            % (+,-)
            L(gfnodes2,gfnodes1) = L(gfnodes2,gfnodes1) - kp*MM2;
            % (-,+)
            L(gfnodes1,gfnodes2) = L(gfnodes1,gfnodes2) - kp*MM1;
            % Gradient Terms
            % ------------------------------------------------------------------
            % (+,+)
            L(gcnodes2,gcnodes2) = L(gcnodes2,gcnodes2) + 0.5*D(2)*(G2 + G2');
            % (-,-)
            L(gcnodes1,gcnodes1) = L(gcnodes1,gcnodes1) - 0.5*D(1)*(G1 + G1');
            % (+,-)
            L(gcnodes2,gcnodes1) = L(gcnodes2,gcnodes1) + 0.5*(D(1)*CG1' - D(2)*CG2);
            % (-,+)
            L(gcnodes1,gcnodes2) = L(gcnodes1,gcnodes2) - 0.5*(D(2)*CG2' - D(1)*CG1);
        end
    % Boundary Face
    else
        matids = mesh.MatID(fcells(1));
        h = mesh.OrthogonalProjection(f,1);
        M = FE.FaceMassMatrix{f,1};
        F = FE.FaceFunctionMatrix{f,1};
        G = FE.FaceGradientMatrix{f,1};
        G = cell_dot(dim,fnorm,G);
        fcnodes = DoF.FaceCellNodes{f,1};
        cnodes = DoF.ConnectivityArray{fcells(1)};
        % Apply boundary terms
        for g=1:ndat.numberEnergyGroups
            gfnodes = fcnodes + (g-1)*ndg;
            gcnodes =  cnodes + (g-1)*ndg;
            D = ndat.Diffusion.DiffXS(matids,g);
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D, h, fflag);
            switch(ndat.Diffusion.BCFlags(fflag))
                case(glob.Dirichlet)
                    L(gfnodes,gfnodes) = L(gfnodes,gfnodes) + kp*M;
                    L(gcnodes,gcnodes) = L(gcnodes,gcnodes) - D*(G + G');
                    rhs(gfnodes) = rhs(gfnodes) +  ndat.Diffusion.BCVals(fflag,g)*(kp*M)*ones(length(gfnodes),1);
                    rhs(gcnodes) = rhs(gcnodes) -  ndat.Diffusion.BCVals(fflag,g)*(D*G)*ones(length(gcnodes),1);
                case(glob.Neumann)
                    rhs(gfnodes) = rhs(gfnodes) - ndat.Diffusion.BCVals(fflag,g)*F;
                case(glob.Robin)
                    L(gfnodes,gfnodes) = L(gfnodes,gfnodes) + 0.5*M;
                    rhs(gfnodes) = rhs(gfnodes) + 2*ndat.Diffusion.BCVals(fflag,g)*F;
                case(glob.Function)
                    % matrix contribution
                    L(gfnodes,gfnodes) = L(gfnodes,gfnodes) + kp*M;
                    L(gcnodes,gcnodes) = L(gcnodes,gcnodes) - D*(G + G');
                    % rhs contribution
                    qx = FE.FaceQuadNodes{f,1};
                    qw = FE.FaceQuadWeights{f,1};
                    cb = FE.FaceBasisValues{f,1}';
                    gb = FE.FaceBasisGrads{f,1};
                    ex_sol = ndat.Diffusion.BCVals{fflag,g}(qx);
                    qwx = qw.*ex_sol;
                    rhs(gcnodes) = rhs(gcnodes) + kp*cb*qwx;
                    for q=1:length(qw)
                        rhs(gcnodes) = rhs(gcnodes) - D*qwx(q)*gb(:,:,q)*tfnorm;
                    end
            end
        end
    end
end
L = sparse(L);
% if glob.print_info, fprintf(rev_str); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,rhs] = get_sparse_matrices(x, ndat, mesh, DoF, FE)
global glob
dim = mesh.Dimension;
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
C_IP = ndat.IP_Constant;
% Allocate Memory
rhs = zeros(ndof,1);
I = [];
J = [];
TMAT = [];
rev_str = [];
% Loop through cells
for tcell=1:mesh.TotalCells
    % Print Current Cell Information
%     if glob.print_info
%         msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    onesnodes = ones(ncnodes,1);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        rows = onesnodes*sg;
        cols = (onesnodes*sg)';
        tmat = (ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M);
        I = [I;rows(:)]; J = [J;cols(:)]; TMAT = [TMAT;tmat(:)];
        if ndat.Diffusion.MMS
            gfunc = ndat.Diffusion.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx));
        else
            tvec = ndat.Diffusion.ExtSource(matID,g)*M*ones(ncnodes,1);
        end
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (gg-1)*ndg;
            xsec = ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg) + ndat.Diffusion.ScatteringXS(matID,gg,g);
            tvec = tvec + xsec*M*x(sgg);
        end
    end
    rhs(sg) = rhs(sg) + tvec;
end

% Loop through faces
for f=1:mesh.TotalFaces
    % Print Current Face Information
%     if glob.print_info
%         msg = sprintf('   -> Building Face: %d of %d',f,mesh.TotalFaces);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    fflag = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:); tfnorm = fnorm';
    % Interior Face
    if fflag == 0
        matids = mesh.MatID(fcells);
        D = ndat.Diffusion.DiffXS(matids,:)';
        h = mesh.OrthogonalProjection(f,:);
        fcnodes1 = DoF.FaceCellNodes{f,1}; fcnodes2 = DoF.FaceCellNodes{f,2};
        cnodes1 = DoF.ConnectivityArray{fcells(1)}; cnodes2 = DoF.ConnectivityArray{fcells(2)};
        conesnodes1 = ones(length(cnodes1),1); conesnodes2 = ones(length(cnodes2),1);
        M1 = FE.FaceMassMatrix{f,1}; M2 = FE.FaceMassMatrix{f,2};
        MM1 = FE.FaceConformingMassMatrix{f,1}; MM2 = FE.FaceConformingMassMatrix{f,2};
        G1 = FE.FaceGradientMatrix{f,1}; G1 = cell_dot(dim,fnorm,G1);
        G2 = FE.FaceGradientMatrix{f,2}; G2 = cell_dot(dim,fnorm,G2);
        CG1 = FE.FaceCouplingGradientMatrix{f,1}; CG1 = cell_dot(dim,fnorm,CG1);
        CG2 = FE.FaceCouplingGradientMatrix{f,2}; CG2 = cell_dot(dim,fnorm,CG2);
        fonesnodes = ones(length(fcnodes1),1);
        % Apply Interior Terms
        for g=1:ndat.numberEnergyGroups
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D(g,:), h, fflag);
            gfnodes1 = fcnodes1 + (g-1)*ndg;
            gfnodes2 = fcnodes2 + (g-1)*ndg;
            gcnodes1 =  cnodes1 + (g-1)*ndg;
            gcnodes2 =  cnodes2 + (g-1)*ndg;
            % Cell rows/columns
            crows11 = conesnodes1*gcnodes1; ccols11 = (conesnodes1*gcnodes1)';
            crows22 = conesnodes2*gcnodes2; ccols22 = (conesnodes2*gcnodes2)';
            crows12 = conesnodes1*gcnodes2; ccols12 = (conesnodes1*gcnodes2)';
            crows21 = conesnodes2*gcnodes1; ccols21 = (conesnodes2*gcnodes1)';
            % Face rows/columns
            frows1 = fonesnodes*gfnodes1; fcols1 = (fonesnodes*gfnodes1)';
            frows2 = fonesnodes*gfnodes2; fcols2 = (fonesnodes*gfnodes2)';
            % Mass Terms
            % ------------------------------------------------------------------
            % (-,-)
            I = [I;frows1(:)]; J = [J;fcols1(:)];
            tmat = kp*M1; TMAT = [TMAT;tmat(:)];
            % (+,+)
            I = [I;frows2(:)]; J = [J;fcols2(:)];
            tmat = kp*M2; TMAT = [TMAT;tmat(:)];
            % (+,-)
            I = [I;frows1(:)]; J = [J;fcols2(:)];
            tmat = -kp*MM2; TMAT = [TMAT;tmat(:)];
            % (-,+)
            I = [I;frows2(:)]; J = [J;fcols1(:)];
            tmat = -kp*MM1; TMAT = [TMAT;tmat(:)];
            % Gradient Terms
            % ------------------------------------------------------------------
            % (-,-)
            I = [I;crows11(:)]; J = [J;ccols11(:)];
            tmat = -0.5*D(1)*(G1 + G1'); TMAT = [TMAT;tmat(:)];
            % (+,+)
            I = [I;crows22(:)]; J = [J;ccols22(:)];
            tmat =  0.5*D(2)*(G2 + G2'); TMAT = [TMAT;tmat(:)];
            % (-,+)
            I = [I;crows21(:)]; J = [J;ccols12(:)];
            tmat =  0.5*(D(1)*CG1' - D(2)*CG2); TMAT = [TMAT;tmat(:)];
            % (+,-)
            I = [I;crows12(:)]; J = [J;ccols21(:)];
            tmat = -0.5*(D(2)*CG2' - D(1)*CG1); TMAT = [TMAT;tmat(:)];
        end
    % Boundary Face
    else
        matids = mesh.MatID(fcells(1));
        h = mesh.OrthogonalProjection(f,1);
        M = FE.FaceMassMatrix{f,1};
        G = FE.FaceGradientMatrix{f,1};
        G = cell_dot(dim,fnorm,G);
        fcnodes = DoF.FaceCellNodes{f,1};
        cnodes = DoF.ConnectivityArray{fcells(1)};
        fonesnodes = ones(length(fcnodes),1);
        conesnodes = ones(length(cnodes),1);
        % Apply boundary terms
        for g=1:ndat.numberEnergyGroups
            gfnodes = fcnodes + (g-1)*ndg;
            gcnodes =  cnodes + (g-1)*ndg;
            crows = conesnodes*gcnodes; ccols = (conesnodes*gcnodes)';
            frows = fonesnodes*gfnodes; fcols = (fonesnodes*gfnodes)';
            D = ndat.Diffusion.DiffXS(matids,g);
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D, h, fflag);
            switch(ndat.Diffusion.BCFlags(fflag))
                case(glob.Dirichlet)
                    tcmat = -D*(G + G'); tfmat = kp*M;
                    I = [I;crows(:)]; J = [J;ccols(:)]; TMAT = [TMAT;tcmat(:)];
                    I = [I;frows(:)]; J = [J;fcols(:)]; TMAT = [TMAT;tfmat(:)];
                    rhs(gfnodes) = rhs(gfnodes) +  ndat.Diffusion.BCVals(fflag,g)*(kp*M)*ones(length(gfnodes),1);
                    rhs(gcnodes) = rhs(gcnodes) -  ndat.Diffusion.BCVals(fflag,g)*(D*G)*ones(length(gcnodes),1);
                case(glob.Neumann)
                    rhs(gfnodes) = rhs(gfnodes) - ndat.Diffusion.BCVals(fflag,g)*M*ones(length(gfnodes),1);
                case(glob.Robin)
                    tmat = 0.5*M;
                    I = [I;frows(:)]; J = [J;fcols(:)]; TMAT = [TMAT;tmat(:)];
                    rhs(gfnodes) = rhs(gfnodes) + 2*ndat.Diffusion.BCVals(fflag,g)*M*ones(length(gfnodes),1);
                case(glob.Function)
                    tcmat = -D*(G + G'); tfmat = kp*M;
                    I = [I;crows(:)]; J = [J;ccols(:)]; TMAT = [TMAT;tcmat(:)];
                    I = [I;frows(:)]; J = [J;fcols(:)]; TMAT = [TMAT;tfmat(:)];
                    qx = FE.FaceQuadNodes{f,1};
                    qw = FE.FaceQuadWeights{f,1};
                    cb = FE.FaceBasisValues{f,1}';
                    gb = FE.FaceBasisGrads{f,1};
                    ex_sol = ndat.Diffusion.BCVals{fflag,g}(qx);
                    qwx = qw.*ex_sol;
                    rhs(gcnodes) = rhs(gcnodes) + kp*cb*qwx;
                    for q=1:length(qw)
                        rhs(gcnodes) = rhs(gcnodes) - D*qwx(q)*gb(:,:,q)*tfnorm;
                    end
            end
        end
    end
end
L = sparse(I,J,TMAT,ndof,ndof);
% if glob.print_info, fprintf(rev_str); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = get_rhs(x, ndat, mesh, DoF, FE)
global glob
dim = mesh.Dimension;
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
C_IP = ndat.IP_Constant;
% Allocate Memory
rhs = zeros(ndof,1);
rev_str = [];
% Loop through cells
for tcell=1:mesh.TotalCells
    % Print Current Cell Information
%     if glob.print_info
%         msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end 
    
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        if ndat.Diffusion.MMS
            gfunc = ndat.Diffusion.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx));
        else
            tvec = ndat.Diffusion.ExtSource(matID,g)*M*ones(ncnodes,1);
        end
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (gg-1)*ndg;
            xsec = ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg) + ndat.Diffusion.ScatteringXS(matID,gg,g);
            tvec = tvec + xsec*M*x(sgg);
        end
    end
    rhs(sg) = rhs(sg) + tvec;
end

% Loop through faces
for f=1:mesh.TotalFaces
    % Print Current Face Information
%     if glob.print_info
%         msg = sprintf('   -> Building Face: %d of %d',f,mesh.TotalFaces);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    fflag = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:); tfnorm = fnorm';
    % Boundary Face
    if fflag ~= 0
        matids = mesh.MatID(fcells(1));
        h = mesh.OrthogonalProjection(f,1);
        M = FE.FaceMassMatrix{f,1};
        G = FE.FaceGradientMatrix{f,1};
        G = cell_dot(dim,fnorm,G);
        fcnodes = DoF.FaceCellNodes{f,1};
        cnodes = DoF.ConnectivityArray{fcells(1)};
        % Apply boundary terms
        for g=1:ndat.numberEnergyGroups
            gfnodes = fcnodes + (g-1)*ndg;
            gcnodes =  cnodes + (g-1)*ndg;
            D = ndat.Diffusion.DiffXS(matids,g);
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D, h, fflag);
            switch(ndat.Diffusion.BCFlags(fflag))
                case(glob.Dirichlet)
                    rhs(gfnodes) = rhs(gfnodes) +  ndat.Diffusion.BCVals(fflag,g)*(kp*M)*ones(length(gfnodes),1);
                    rhs(gcnodes) = rhs(gcnodes) -  ndat.Diffusion.BCVals(fflag,g)*(D*G)*ones(length(gcnodes),1);
                case(glob.Neumann)
                    rhs(gfnodes) = rhs(gfnodes) - ndat.Diffusion.BCVals(fflag,g)*M*ones(length(gfnodes),1);
                case(glob.Robin)
                    rhs(gfnodes) = rhs(gfnodes) + 2*ndat.Diffusion.BCVals(fflag,g)*M*ones(length(gfnodes),1);
                case(glob.Function)
                    qx = FE.FaceQuadNodes{f,1};
                    qw = FE.FaceQuadWeights{f,1};
                    cb = FE.FaceBasisValues{f,1}';
                    gb = FE.FaceBasisGrads{f,1};
                    ex_sol = ndat.Diffusion.BCVals{fflag,g}(qx);
                    qwx = qw.*ex_sol;
                    rhs(gcnodes) = rhs(gcnodes) + kp*cb*qwx;
                    for q=1:length(qw)
                        rhs(gcnodes) = rhs(gcnodes) - D*qwx(q)*gb(:,:,q)*tfnorm;
                    end
            end
        end
    end
end
% if glob.print_info, fprintf(rev_str); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_Ax(x, ndat, mesh, DoF, FE)
global glob
dim = mesh.Dimension;
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
C_IP = ndat.IP_Constant;
% Allocate Memory
out = zeros(ndof,1);
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        out(sg) = out(sg) + (ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M)*x(sg);
    end
end
% Loop through faces
for f=1:mesh.TotalFaces
    fflag = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:);
    % Interior Face
    if fflag == 0
        matids = mesh.MatID(fcells);
        D = ndat.Diffusion.DiffXS(matids,:)';
        h = mesh.OrthogonalProjection(f,:);
        fcnodes1 = DoF.FaceCellNodes{f,1}; fcnodes2 = DoF.FaceCellNodes{f,2};
        cnodes1 = DoF.ConnectivityArray{fcells(1)}; cnodes2 = DoF.ConnectivityArray{fcells(2)};
        M1 = FE.FaceMassMatrix{f,1}; M2 = FE.FaceMassMatrix{f,2};
        MM1 = FE.FaceConformingMassMatrix{f,1}; MM2 = FE.FaceConformingMassMatrix{f,2};
        G1 = FE.FaceGradientMatrix{f,1}; G1 = cell_dot(dim,fnorm,G1);
        G2 = FE.FaceGradientMatrix{f,2}; G2 = cell_dot(dim,fnorm,G2);
        CG1 = FE.FaceCouplingGradientMatrix{f,1}; CG1 = cell_dot(dim,fnorm,CG1);
        CG2 = FE.FaceCouplingGradientMatrix{f,2}; CG2 = cell_dot(dim,fnorm,CG2);
        % Apply Interior Terms
        for g=1:ndat.numberEnergyGroups
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D(g,:), h, fflag);
            gfnodes1 = fcnodes1 + (g-1)*ndg;
            gfnodes2 = fcnodes2 + (g-1)*ndg;
            gcnodes1 =  cnodes1 + (g-1)*ndg;
            gcnodes2 =  cnodes2 + (g-1)*ndg;
            % Mass Terms
            % ------------------------------------------------------------------
            % (-,-)
            out(gfnodes1) = out(gfnodes1) + kp*M1*x(gfnodes1);
            % (+,+)
            out(gfnodes2) = out(gfnodes2) + kp*M2*x(gfnodes2);
            % (+,-)
            out(gfnodes2) = out(gfnodes2) - kp*MM2*x(gfnodes1);
            % (-,+)
            out(gfnodes1) = out(gfnodes1) - kp*MM1*x(gfnodes2);
            % Gradient Terms
            % ------------------------------------------------------------------
            % (+,+)
            out(gcnodes2) = out(gcnodes2) + 0.5*D(2)*(G2 + G2')*x(gcnodes2);
            % (-,-)
            out(gcnodes1) = out(gcnodes1) - 0.5*D(1)*(G1 + G1')*x(gcnodes1);
            % (+,-)
            out(gcnodes2) = out(gcnodes2) + 0.5*(D(1)*CG1' - D(2)*CG2)*x(gcnodes1);
            % (-,+)
            out(gcnodes1) = out(gcnodes1) - 0.5*(D(2)*CG2' - D(1)*CG1)*x(gcnodes2);
        end
    % Boundary Face
    else
        matids = mesh.MatID(fcells(1));
        h = mesh.OrthogonalProjection(f,1);
        M = FE.FaceMassMatrix{f,1};
        G = FE.FaceGradientMatrix{f,1};
        G = cell_dot(dim,fnorm,G);
        fcnodes = DoF.FaceCellNodes{f,1};
        cnodes = DoF.ConnectivityArray{fcells(1)};
        % Apply boundary terms
        for g=1:ndat.numberEnergyGroups
            gfnodes = fcnodes + (g-1)*ndg;
            gcnodes =  cnodes + (g-1)*ndg;
            D = ndat.Diffusion.DiffXS(matids,g);
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D, h, fflag);
            switch(ndat.Diffusion.BCFlags(fflag))
                case(glob.Dirichlet)
                    out(gfnodes) = out(gfnodes) + kp*M*x(gfnodes);
                    out(gcnodes) = out(gcnodes) - D*(G + G')*x(gcnodes);
                case(glob.Robin)
                    out(gfnodes) = out(gfnodes) + 0.5*M*x(gfnodes);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_penalty_coefficient(C,p,D,h,eflag)
c = C*(1+p)*p;
if eflag == 0
    out = c/2*(D(1)/h(1) + D(2)/h(2));
else
    out = c*D/h;
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