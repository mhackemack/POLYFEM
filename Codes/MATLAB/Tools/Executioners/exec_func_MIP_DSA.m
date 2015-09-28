%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get MIP DSA System Matrix
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Retrieves the MIP system matrix for a 1-group problem. It is
%                   assumed that the cross sections have been properly handled
%                   prior to this function call.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = exec_func_MIP_DSA(data,accel_id,xsid,mesh,DoF,FE)
global glob
% Retrieve MIP DSA System Matrix
if DoF.TotalDoFs < glob.maxMatrix
    A = get_global_matrices(data.XS(xsid),data.Acceleration.Info(accel_id),mesh,DoF,FE);
else
    A = get_sparse_matrices(data.XS(xsid),data.Acceleration.Info(accel_id),mesh,DoF,FE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = get_global_matrices(XS, Accel, mesh, DoF, FE)
global glob
dim = mesh.Dimension;
C_IP = Accel.IP_Constant;
% Allocate Memory
L = zeros(DoF.TotalDoFs,DoF.TotalDoFs);
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    L(cnodes,cnodes) = L(cnodes,cnodes) + XS.DiffXS(matID)*K + XS.AbsorbXS(matID)*M;
end
% Loop through faces
for f=1:mesh.TotalFaces
    fflag = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:);
    % Interior Face
    if fflag == 0
        matids = mesh.MatID(fcells);
        D = XS.DiffXS(matids);
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
        kp = get_penalty_coefficient(C_IP, DoF.Degree, D, h, fflag);
        % Mass Terms
        % -------------------------------------------------------------
        % (-,-)
        L(fcnodes1,fcnodes1) = L(fcnodes1,fcnodes1) + kp*M1;
        % (+,+)
        L(fcnodes2,fcnodes2) = L(fcnodes2,fcnodes2) + kp*M2;
        % (+,-)
        L(fcnodes2,fcnodes1) = L(fcnodes2,fcnodes1) - kp*MM2;
        % (-,+)
        L(fcnodes1,fcnodes2) = L(fcnodes1,fcnodes2) - kp*MM1;
        % Gradient Terms
        % -------------------------------------------------------------
        % (+,+)
        L(cnodes2,cnodes2) = L(cnodes2,cnodes2) + 0.5*D(2)*(G2 + G2');
        % (-,-)
        L(cnodes1,cnodes1) = L(cnodes1,cnodes1) - 0.5*D(1)*(G1 + G1');
        % (+,-)
        L(cnodes2,cnodes1) = L(cnodes2,cnodes1) + 0.5*(D(1)*CG1' - D(2)*CG2);
        % (-,+)
        L(cnodes1,cnodes2) = L(cnodes1,cnodes2) - 0.5*(D(2)*CG2' - D(1)*CG1);
    % Boundary Face
    else
        matids = mesh.MatID(fcells(1));
        D = XS.DiffXS(matids);
        h = mesh.OrthogonalProjection(f,1);
        M = FE.FaceMassMatrix{f,1};
        G = FE.FaceGradientMatrix{f,1};
        G = cell_dot(dim,fnorm,G);
        fcnodes = DoF.FaceCellNodes{f,1};
        cnodes = DoF.ConnectivityArray{fcells(1)};
        kp = get_penalty_coefficient(C_IP, DoF.Degree, D, h, fflag);
        if     (XS.BCFlags(fflag) == glob.Vacuum || ...
                XS.BCFlags(fflag) == glob.IncidentIsotropic || ...
                XS.BCFlags(fflag) == glob.IncidentCurrent || ...
                XS.BCFlags(fflag) == glob.IncidentBeam)
            L(fcnodes,fcnodes) = L(fcnodes,fcnodes) + kp*M;
            L(cnodes,cnodes) = L(cnodes,cnodes) - 0.5*D*(G + G');
        end
    end
end
L = sparse(L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,rhs] = get_sparse_matrices(XS, Accel, mesh, DoF, FE)
global glob
dim = mesh.Dimension;
C_IP = Accel.IP_Constant;
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
    rows = onesnodes*cnodes;
    cols = (onesnodes*cnodes)';
    tmat = XS.DiffXS(matID)*K + XS.AbsorbXS(matID)*M;
    I = [I;rows(:)]; J = [J;cols(:)]; TMAT = [TMAT;tmat(:)];
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
        D = zeros(2,1);
        h = zeros(2,1);
        matids = mesh.MatID(fcells);
        for c=1:2
            h(c) = mesh.OrthogonalProjection(f,c);
            fcnodes{c} = DoF.FaceCellNodes{f,c};
            cnodes{c} = DoF.ConnectivityArray{fcells(c)};
            conesnodes{c} = ones(length(cnodes{c}),1);
            D(c) = XS.DiffXS(matids(c));
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
        D = XS.DiffXS(matids);
        h = mesh.OrthogonalProjection(f,1);
        kp = get_penalty_coefficient(C_IP, DoF.Degree, D, h, fflag);
        M = FE.FaceMassMatrix{f,1};
        G = FE.FaceGradientMatrix{f,1};
        G = cell_dot(dim,fnorm,G);
        fcnodes = DoF.FaceCellNodes{f,1};
        cnodes = DoF.ConnectivityArray{fcells(1)};
        fonesnodes = ones(length(fcnodes),1);
        conesnodes = ones(length(cnodes),1);
        % Apply boundary terms
        crows = conesnodes*cnodes; ccols = (conesnodes*cnodes)';
        frows = fonesnodes*fcnodes; fcols = (fonesnodes*fcnodes)';
        if     (XS.BCFlags(fflag) == glob.Vacuum || ...
                XS.BCFlags(fflag) == glob.IncidentIsotropic || ...
                XS.BCFlags(fflag) == glob.IncidentCurrent || ...
                XS.BCFlags(fflag) == glob.IncidentBeam)
            tcmat = -0.5*D*(G + G'); tfmat = kp*M;
            I = [I;crows(:)]; J = [J;ccols(:)]; TMAT = [TMAT;tcmat(:)];
            I = [I;frows(:)]; J = [J;fcols(:)]; TMAT = [TMAT;tfmat(:)];
        end
    end
end
L = sparse(I,J,TMAT,ndof,ndof);
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