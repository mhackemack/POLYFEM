%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve DCF DSA Diffusion Problem
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_DCF_DSA(ndat, solvdat, mesh, DoF, FE, x, A)
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
    rT = solvdat.relativeTolerance;
    mI = solvdat.maxIterations;
    [x,flag,res,it] = pcg(@(x) get_Ax(x, ndat, mesh, DoF, FE),rhs,rT,mI,[],[],x);
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
%                              Function Listing
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
% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------
for f=1:mesh.TotalFaces
    fflag = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:);
    % Interior Face
    if fflag == 0
        matids = mesh.MatID(fcells);
        D = ndat.Diffusion.DiffXS(matids,1);
        fcnodes1 = DoF.FaceCellNodes{f,1}; fcnodes2 = DoF.FaceCellNodes{f,2};
        cnodes1 = DoF.ConnectivityArray{fcells(1)}; cnodes2 = DoF.ConnectivityArray{fcells(2)};
        M1 = FE.FaceMassMatrix{f,1}; M2 = FE.FaceMassMatrix{f,2};
        MM1 = FE.FaceConformingMassMatrix{f,1}; MM2 = FE.FaceConformingMassMatrix{f,2};
        G1 = FE.FaceGradientMatrix{f,1}; G1 = cell_dot(dim,fnorm,G1);
        G2 = FE.FaceGradientMatrix{f,2}; G2 = cell_dot(dim,fnorm,G2);
        CG1 = FE.FaceCouplingGradientMatrix{f,1}; CG1 = cell_dot(dim,fnorm,CG1);
        CG2 = FE.FaceCouplingGradientMatrix{f,2}; CG2 = cell_dot(dim,fnorm,CG2);
        % MIP - Mass Terms
        % ----------------------------------------------------------------------
        L(fcnodes1,fcnodes1) = L(fcnodes1,fcnodes1) + 0.25*M1;
        L(fcnodes2,fcnodes2) = L(fcnodes2,fcnodes2) + 0.25*M2;
        L(fcnodes2,fcnodes1) = L(fcnodes2,fcnodes1) - 0.25*MM2;
        L(fcnodes1,fcnodes2) = L(fcnodes1,fcnodes2) - 0.25*MM1;
        % MIP - Gradient Terms
        % ----------------------------------------------------------------------
        L(cnodes2,cnodes2) = L(cnodes2,cnodes2) + 0.5*D(2)*(G2 + G2');
        L(cnodes1,cnodes1) = L(cnodes1,cnodes1) - 0.5*D(1)*(G1 + G1');
        L(cnodes2,cnodes1) = L(cnodes2,cnodes1) + 0.5*(D(1)*CG1' - D(2)*CG2);
        L(cnodes1,cnodes2) = L(cnodes1,cnodes2) - 0.5*(D(2)*CG2' - D(1)*CG1);
        % Type 3 Edges - 
        % ----------------------------------------------------------------------
        
        % Type 4 Edges - 
        % ----------------------------------------------------------------------
        
    % Boundary Face
    else
        matids = mesh.MatID(fcells(1));
        M = FE.FaceMassMatrix{f,1};
        G = FE.FaceGradientMatrix{f,1};
        G = cell_dot(dim,fnorm,G);
        G3 = FE.FaceGradNormalMatrix{f,1};
        G4 = FE.FaceStiffnessMatrix{f,1};
        fcnodes = DoF.FaceCellNodes{f,1};
        cnodes = DoF.ConnectivityArray{fcells(1)};
        D = ndat.Diffusion.DiffXS(matids,1);
        if     (ndat.Transport.BCFlags(fflag) == glob.Vacuum || ...
                ndat.Transport.BCFlags(fflag) == glob.IncidentIsotropic || ...
                ndat.Transport.BCFlags(fflag) == glob.IncidentCurrent || ...
                ndat.Transport.BCFlags(fflag) == glob.IncidentBeam)
            % MIP Terms
            % ------------------------------------------------------------------
            L(fcnodes,fcnodes) = L(fcnodes,fcnodes) + 0.25*M;
            L(cnodes,cnodes) = L(cnodes,cnodes) - 0.5*D*(G + G');
            % Type 3/4 Edges
            % ------------------------------------------------------------------
            L(cnodes,cnodes) = L(cnodes,cnodes) - 9/16*D*D*(G3+G4);
            % Q1 Terms
            % ------------------------------------------------------------------
            
        elseif  ndat.Transport.BCFlags(fflag) == glob.Reflecting || ...
                ndat.Transport.BCFlags(fflag) == glob.Periodic
            % J/Y Terms
            % ------------------------------------------------------------------
            Jin = ndat.Transport.OutgoingCurrents{f}(:,g) - ndat.Transport.OutgoingCurrentsOld{f}(:,g);
            rhs(fcnodes) = rhs(fcnodes) + M*Jin;
            % Q1 Terms
            % ------------------------------------------------------------------
            
        end
    end
end
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
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    
end
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
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_Ax(x, ndat, mesh, DoF, FE)
global glob
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
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%