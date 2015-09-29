%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve MIP DSA Diffusion Problem
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_MIP_DSA(ndat, solvdat, mesh, DoF, FE, x, A)
% Throw error if opposing reflecting boundaries are present - this will be
% resolved at a later date...maybe
% ------------------------------------------------------------------------------
if ndat.Transport.HasOpposingReflectingBoundary
    error('Currently cannot support opposing reflecting boundaries.');
end
% ------------------------------------------------------------------------------
global glob
nout = nargout;
ndg = DoF.TotalDoFs;
ndof = ndat.numberEnergyGroups * ndg;
if nargin < 5 || isempty(x)
    x = ones(ndof,1);
else
    x = cell_to_vector(x, DoF);
end
if length(x) > glob.maxSparse
    rhs = get_rhs(x, ndat, mesh, DoF, FE);
    [x,flag,res,it] = pcg(@(x) get_Ax(x, ndat, mesh, DoF, FE),rhs,1e-3,1e5,[],[],x);
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
% Outputs
x = vector_to_cell(x,DoF);
varargout{1} = x;
if nout == 2
    if exist('A', 'var')
        varargout{2} = A; 
    else
        varargout{2} = []; 
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
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        L(sg,sg) = L(sg,sg) + ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.AbsorbXS(matID,g)*M;
        rhs(sg) = rhs(sg) + ndat.Diffusion.ScatteringXS(matID,g,g) * M * x(sg);
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
            % -------------------------------------------------------------
            % (-,-)
            L(gfnodes1,gfnodes1) = L(gfnodes1,gfnodes1) + kp*M1;
            % (+,+)
            L(gfnodes2,gfnodes2) = L(gfnodes2,gfnodes2) + kp*M2;
            % (+,-)
            L(gfnodes2,gfnodes1) = L(gfnodes2,gfnodes1) - kp*MM2;
            % (-,+)
            L(gfnodes1,gfnodes2) = L(gfnodes1,gfnodes2) - kp*MM1;
            % Gradient Terms
            % -------------------------------------------------------------
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
            if     (ndat.Transport.BCFlags(fflag) == glob.Vacuum || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentIsotropic || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentCurrent || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentBeam)
                L(gfnodes,gfnodes) = L(gfnodes,gfnodes) + kp*M;
                L(gcnodes,gcnodes) = L(gcnodes,gcnodes) - 0.5*D*(G + G');
            end
        end
    end
end
L = sparse(L);
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
% Loop through cells
for tcell=1:mesh.TotalCells
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
        tmat = (ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.AbsorbXS(matID,g)*M);
        I = [I;rows(:)]; J = [J;cols(:)]; TMAT = [TMAT;tmat(:)];
        rhs(sg) = rhs(sg) + ndat.Diffusion.ScatteringXS(matID,g,g) * M * x(sg);
    end
end
% Loop through faces
for f=1:mesh.TotalFaces
    fflag = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:);
    % Interior Face
    if fflag == 0
        fcnodes = cell(2,1);
        cnodes = cell(2,1);
        M = cell(2,1);
        MM = cell(2,1);
        G = cell(2,1);
        CG = cell(2,1);
        conesnodes = cell(2,1);
        D = zeros(ndat.numberEnergyGroups,2);
        h = zeros(2,1);
        matids = mesh.MatID(fcells);
        for c=1:2
            h(c) = mesh.OrthogonalProjection(f,c);
            fcnodes{c} = DoF.FaceCellNodes{f,c};
            cnodes{c} = DoF.ConnectivityArray{fcells(c)};
            conesnodes{c} = ones(length(cnodes{c}),1);
            for g=1:ndat.numberEnergyGroups
                D(g,c) = ndat.Diffusion.DiffXS(matids(c),g);
            end
            M{c} = FE.FaceMassMatrix{f,c};
            MM{c} = FE.FaceConformingMassMatrix{f,c};
            G{c} = FE.FaceGradientMatrix{f,c};
            G{c} = cell_dot(dim,fnorm,G{c});
            CG{c} = FE.FaceCouplingGradientMatrix{f,c};
            CG{c} = cell_dot(dim,fnorm,CG{c});
        end
        fonesnodes = ones(length(fcnodes{1}),1);
        % Apply Interior Terms
        for g=1:ndat.numberEnergyGroups
            kp = get_penalty_coefficient(C_IP, DoF.Degree, D(g,:), h, fflag);
            gfnodes1 = fcnodes{1} + (g-1)*ndg;
            gfnodes2 = fcnodes{2} + (g-1)*ndg;
            gcnodes1 =  cnodes{1} + (g-1)*ndg;
            gcnodes2 =  cnodes{2} + (g-1)*ndg;
            % Cell rows/columns
            crows11 = conesnodes{1}*gcnodes1; ccols11 = (conesnodes{1}*gcnodes1)';
            crows22 = conesnodes{2}*gcnodes2; ccols22 = (conesnodes{2}*gcnodes2)';
            crows12 = conesnodes{1}*gcnodes2; ccols12 = (conesnodes{1}*gcnodes2)';
            crows21 = conesnodes{2}*gcnodes1; ccols21 = (conesnodes{2}*gcnodes1)';
            % Face rows/columns
            frows1 = fonesnodes*gfnodes1; fcols1 = (fonesnodes*gfnodes1)';
            frows2 = fonesnodes*gfnodes2; fcols2 = (fonesnodes*gfnodes2)';
            % Mass Terms
            % ------------------------------------------------------------------
            % (-,-)
            I = [I;frows1(:)]; J = [J;fcols1(:)];
            tmat = kp*M{1}; TMAT = [TMAT;tmat(:)];
            % (+,+)
            I = [I;frows2(:)]; J = [J;fcols2(:)];
            tmat = kp*M{2}; TMAT = [TMAT;tmat(:)];
            % (+,-)
            I = [I;frows1(:)]; J = [J;fcols2(:)];
            tmat = -kp*MM{2}; TMAT = [TMAT;tmat(:)];
            % (-,+)
            I = [I;frows2(:)]; J = [J;fcols1(:)];
            tmat = -kp*MM{1}; TMAT = [TMAT;tmat(:)];
            % Gradient Terms
            % ------------------------------------------------------------------
            % (-,-)
            I = [I;crows11(:)]; J = [J;ccols11(:)];
            tmat = -0.5*D(1)*(G{1} + G{1}'); TMAT = [TMAT;tmat(:)];
            % (+,+)
            I = [I;crows22(:)]; J = [J;ccols22(:)];
            tmat =  0.5*D(2)*(G{2} + G{2}'); TMAT = [TMAT;tmat(:)];
            % (-,+)
            I = [I;crows21(:)]; J = [J;ccols12(:)];
            tmat =  0.5*(D(1)*CG{1}' - D(2)*CG{2}); TMAT = [TMAT;tmat(:)];
            % (+,-)
            I = [I;crows12(:)]; J = [J;ccols21(:)];
            tmat = -0.5*(D(2)*CG{2}' - D(1)*CG{1}); TMAT = [TMAT;tmat(:)];
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
            if     (ndat.Transport.BCFlags(fflag) == glob.Vacuum || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentIsotropic || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentCurrent || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentBeam)
                tcmat = -0.5*D*(G + G'); tfmat = kp*M;
                I = [I;crows(:)]; J = [J;ccols(:)]; TMAT = [TMAT;tcmat(:)];
                I = [I;frows(:)]; J = [J;fcols(:)]; TMAT = [TMAT;tfmat(:)];
            end
        end
    end
end
L = sparse(I,J,TMAT,ndof,ndof);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = get_rhs(x, ndat, mesh, DoF, FE)
ndg = DoF.TotalDoFs;
ng = ndat.numberEnergyGroups;
ndof = ng*ndg;
% Allocate Memory
rhs = zeros(ndof,1);
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        rhs(sg) = rhs(sg) + ndat.Diffusion.ScatteringXS(matID,g,g) * M * x(sg);
    end
end
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
        out(sg) = out(sg) + (ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.AbsorbXS(matID,g)*M)*x(sg);
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
            if     (ndat.Transport.BCFlags(fflag) == glob.Vacuum || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentIsotropic || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentCurrent || ...
                    ndat.Transport.BCFlags(fflag) == glob.IncidentBeam)
                out(gfnodes) = out(gfnodes) + kp*M*x(gfnodes);
                out(gcnodes) = out(gcnodes) - 0.5*D*(G + G')*x(gcnodes);
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
out = max(out, 0.25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = cell_dot(dim,vec1, vec2)
if dim == 1
    out = vec1*vec2{1};
elseif dim == 2
    out = vec1(1)*vec2{1} + vec1(2)*vec2{2};
else
    out = vec1(1)*vec2{1} + vec1(2)*vec2{2} + vec1(3)*vec2{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%