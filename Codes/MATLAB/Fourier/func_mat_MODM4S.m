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
function A = func_mat_MODM4S(lam, input)
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
A = zeros(ndofs); %B = zeros(ndofs);
% Loop through Cells and Build Matrices
% -------------------------------------
for c=1:mesh.TotalCells
    cn  = dof.ConnectivityArray{c};
    mat = mesh.MatID(c);
    M   = fe.CellMassMatrix{c};
    K   = fe.CellStiffnessMatrix{c};
    axs = data.AveAbsorbXS(mat);
    D   = data.AveDiffusionXS(mat);
    A(cn,cn) = A(cn,cn) + axs*M + D*K;
end
% Apply Volumetric Phase Shift
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
        mID   = mesh.MatID(fcells);
        D     = data.AveDiffusionXS(mID);
        h     = mesh.OrthogonalProjection(f,:);
        k     = get_penalty_coefficient(data.IP_Constant,fe.Degree,0,D,h);
%         k     = 1/4;
        fn1   = dof.FaceCellNodes{f,1};
        fn2   = dof.FaceCellNodes{f,2};
        cn1   = dof.ConnectivityArray{fcells(1)};
        cn2   = dof.ConnectivityArray{fcells(2)};
        M1    = fe.FaceMassMatrix{f,1}; 
        M2    = fe.FaceMassMatrix{f,2};
        MM1   = fe.FaceConformingMassMatrix{f,1};
        MM2   = fe.FaceConformingMassMatrix{f,2};
        G1    = cell_dot(dim,fnorm,fe.FaceGradientMatrix{f,1});
        G2    = cell_dot(dim,fnorm,fe.FaceGradientMatrix{f,2});
        CG1    = cell_dot(dim,fnorm,fe.FaceCouplingGradientMatrix{f,1});
        CG2    = cell_dot(dim,fnorm,fe.FaceCouplingGradientMatrix{f,2});
        PMf1  = PM(fn1, fn1); PMf2 = PM(fn2, fn2);
        PMc1  = PM(cn1, cn1); PMc2 = PM(cn2, cn2);
        % Build all interior face-coupling matrix contributions
        % -----------------------------------------------------
        % ( [[b]] , [[u]] )
        A(fn1,fn1) = A(fn1,fn1) + k*M1*PMf1;   % (-,-)
        A(fn1,fn2) = A(fn1,fn2) - k*MM1*PMf2;  % (-,+)
        A(fn2,fn1) = A(fn2,fn1) - k*MM2*PMf1;  % (+,-)
        A(fn2,fn2) = A(fn2,fn2) + k*M2*PMf2;   % (+,+)
        % ( [[b]] , {{Du}} )
        A(cn1,cn1) = A(cn1,cn1) - D(1)/2*G1'*PMc1;   % (-,-)
        A(cn1,cn2) = A(cn1,cn2) - D(2)/2*CG2'*PMc2;  % (-,+)
        A(cn2,cn1) = A(cn2,cn1) + D(1)/2*CG1'*PMc1;  % (+,-)
        A(cn2,cn2) = A(cn2,cn2) + D(2)/2*G2'*PMc2;   % (+,+)
    % Boundary Faces
    else
        % Get preliminaries
        op_f  = mesh.PeriodicOppositeFaces(f);
        op_c  = mesh.FaceCells(op_f,1);
        mID   = mesh.MatID([fcells(1), op_c]);
        D     = data.AveDiffusionXS(mID);
        h     = mesh.OrthogonalProjection([f,op_f],1);
%         k     = 1/4;
        k     = get_penalty_coefficient(data.IP_Constant,fe.Degree,0,D,h);
        fn1   = dof.FaceCellNodes{f,1};
        fn2   = dof.PeriodicFaceDoFs{f}(:,2)';
        cn1   = dof.ConnectivityArray{fcells(1)};
        cn2   = dof.ConnectivityArray{op_c};
        M     = fe.FaceMassMatrix{f,1};
        fpnorm= (-1)*fnorm;
        G1    = cell_dot(dim,fnorm,fe.FaceGradientMatrix{f,1});
        G2    = cell_dot(dim,fnorm,fe.FaceGradientMatrix{op_f,1});
        GG1   = cell_dot(dim,fpnorm,fe.FaceGradientMatrix{f,1});
        GG2   = cell_dot(dim,fpnorm,fe.FaceGradientMatrix{op_f,1});
        t_off = zeros(1,dim); t_off(1,off(f,1)) = off(f,2);
        PMf1  = PM(fn1, fn1); PMf2 = PM(fn2, fn2)*exp(1i*t_off*lam);
        PMc1  = PM(cn1, cn1); PMc2 = PM(cn2, cn2)*exp(1i*t_off*lam);
        % Build penalty contributions
        % ----------------------------------------------------------------------
        % ( [[b]] , [[u]] )
            A(fn1,fn1) = A(fn1,fn1) + k*M*PMf1;
            A(fn1,fn2) = A(fn1,fn2) - k*M*PMf2;
        % Build (-) side contributions
        % ----------------------------------------------------------------------
        if sum(fnorm) > 0
            % ( [[b]] , {{Du}} )
            A(cn1,cn1) = A(cn1,cn1) - 0.5*D(1)*G1'*PMc1;
            A(cn1,cn2) = A(cn1,cn2) - 0.5*D(2)*G1'*PMc2;
        % Build (+) side contributions
        % ----------------------------------------------------------------------
        else
            % ( [[b]] , {{Du}} )
            A(cn1,cn1) = A(cn1,cn1) + 0.5*D(1)*GG1'*PMc1;
            A(cn1,cn2) = A(cn1,cn2) + 0.5*D(2)*GG1'*PMc2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxiallary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_penalty_coefficient(c, p, eid , D, h)
C = c*p*(p+1);
if eid == 0
    out = C/2*(D(1)/h(1) + D(2)/h(2));
else
    out = C*D(1)/h(1);
end
out = max(out,0.25);
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
