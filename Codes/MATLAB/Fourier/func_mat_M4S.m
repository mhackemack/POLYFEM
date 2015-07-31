%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate MIP Matrices
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = func_mat_M4S(lam, input)
global glob
% Copy Input Space
% ----------------
data = input.data;
mesh = input.mesh;
dof  = input.dof;
fe   = input.fe;
off  = input.offset;
% Retrieve Preliminary Data
% -------------------------
dim = mesh.Dimension;
ndofs = dof.TotalDoFs;
node_locs = dof.NodeLocations;
if dim == size(lam,2); lam=lam'; end
PV = exp(1i*node_locs*lam);
PM = diag(PV);
% Allocate Matrix Arrays
% ----------------------
A = zeros(ndofs); B = zeros(ndofs);
% Loop through Cells and Build Matrices
% -------------------------------------
for c=1:mesh.TotalCells
    cn  = dof.ConnectivityArray{c};
    mat = mesh.MatID(c);
    M   = fe.CellMassMatrix{c};
    K   = fe.CellStiffnessMatrix{c};
    sxs = data.ScatteringXS(mat);
    axs = data.AbsorbXS(mat);
    D   = data.DiffusionXS(mat);
    B(cn,cn) = B(cn,cn) + sxs*M;
    A(cn,cn) = A(cn,cn) + axs*M + D*K;
end
% Apply Volumetric Phase Shift
B = B * PM;
A = A * PM;
% Loop through Faces and Build Matrices
% -------------------------------------
for f=1:mesh.TotalFaces
    fid = mesh.FaceID(f);
    fcells = mesh.FaceCells(f,:);
    fnorm = mesh.FaceNormal(f,:);
    % Interior Faces
    if fid == 0
        % Get preliminaries
        mID = mesh.MatID(fcells);
        D   = data.DiffusionXS(mID);
        h   = mesh.OrthogonalProjection(f,:);
        fn1 = dof.FaceCellNodes{f,1};
        fn2 = dof.FaceCellNodes{f,2};
        cn1 = dof.ConnectivityArray{fcells(1)};
        cn2 = dof.ConnectivityArray{fcells(2)};
        M   = fe.FaceMassMatrix{f,1};
        MM  = fe.FaceConformingMassMatrix{f,1};
        G1  = cell_dot(dim,fnorm,fe.FaceGradientMatrix{f,1});
        G2  = cell_dot(dim,fnorm,fe.FaceGradientMatrix{f,2});
        PMf1 = PM(fn1, fn1); PMf2 = PM(fn2, fn2);
        PMc1 = PM(cn1, cn1); PMc2 = PM(cn2, cn2);
        % Build all interior face-coupling matrix contributions
        % -----------------------------------------------------
        % ( [[u]] , [[b]] )
        
        % ( {{Du}} , [[b]] )
        
        % ( [[u]] , {{Db}} )
        
    % Boundary Faces
    else
        % Get preliminaries
        op_f  = mesh.PeriodicOppositeFaces(f);
        op_c  = mesh.FaceCells(op_f,1);
        mID   = mesh.MatID([fcells(1), op_c]);
        D     = data.DiffusionXS(mID);
        h     = mesh.OrthogonalProjection([f,op_f],1);
        fn1   = dof.FaceCellNodes{f,1};
        fn2   = dof.PeriodicFaceDoFs{f}(:,2)';
        cn1   = dof.ConnectivityArray{fcells(1)};
        cn2   = dof.ConnectivityArray{op_c};
        M     = fe.FaceMassMatrix{f,1};
        MM    = fe.FaceConformingMassMatrix{f,1};
        G1    = cell_dot(dim,fnorm,fe.FaceGradientMatrix{f,1});
        G2    = cell_dot(dim,fnorm,fe.FaceGradientMatrix{op_f,1});
        t_off = zeros(1:dim); t_off(1,off(f,1)) = off(f,2);
        PMf1  = PM(fn1, fn1); PMf2 = PM(fn2, fn2)*exp(1i*t_off*lam);
        PMc1  = PM(cn1, cn1); PMc2 = PM(cn2, cn2)*exp(1i*t_off*lam);
        % Build all periodic matrix contributions
        % ---------------------------------------
        % ( [[u]] , [[b]] )
        A(fn1,fn1) = A(fn1,fn1) + 0.25*M*PMf1;
        A(fn1,fn2) = A(fn1,fn2) - 0.25*M*PMf2;
        % ( {{Du}} , [[b]] )
        
        % ( [[u]] , {{Db}} )
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxiallary Functions
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
