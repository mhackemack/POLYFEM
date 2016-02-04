%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate SI Matrix
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Generates the flux moment LHS system matrix, T.
%
%                   T = D*L^(-1)*M*S
%
%                   D = Discrete-to-Moment Operator
%                   L = Transport Operator (streaming + interaction)
%                   M = Moment-to-Discrete Operator
%                   S = Scattering Operator
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = func_mat_SI_upwind(lam, input)
% Copy Input Space
% ------------------------------------------------------------------------------
data = input.data;
mesh = input.mesh;
dof = input.dof;
fe = input.fe;
m_quad = input.Quadrature;
m2d = m_quad.moment_to_discrete;
d2m = m_quad.discrete_to_moment;
lamt = lam';
% Retrieve Preliminary Data
% ------------------------------------------------------------------------------
dim = mesh.Dimension;
ndofs = dof.TotalDoFs;
zn = zeros(ndofs);
node_locs = dof.NodeLocations;
if dim == size(lam,2); lam=lam'; end
num_dirs = m_quad.NumberAngularDirections;
PV = exp(1i*node_locs*lam);
PM = diag(PV);
% Allocate Matrix Arrays
% ------------------------------------------------------------------------------
L = cell(num_dirs, 1);
S = zeros(ndofs); T = zeros(ndofs);
for q=1:num_dirs
    L{q} = zeros(ndofs);
end
% Loop through Cells and Build Volumetric Portions of Matrices
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cn  = dof.ConnectivityArray{c};
    mat = mesh.MatID(c);
    M   = fe.CellMassMatrix{c};
    G   = fe.CellGradientMatrix{c};
    txs = data.TotalXS(mat);
    sxs = data.ScatteringXS(mat);
    S(cn,cn) = S(cn,cn) + sxs*M;
    % Loop through quadrature
    for q=1:num_dirs
        GG = cell_dot(dim, m_quad.AngularDirections(q,:), G)';
        L{q}(cn,cn) = L{q}(cn,cn) + txs*M - GG;
    end
end
% Apply Volumetric Phase Shift
% ------------------------------------------------------------------------------
S = S * PM;
for q=1:num_dirs
    L{q} = L{q} * PM;
end
% Build Face Contributions to Transport Matrix
% ------------------------------------------------------------------------------
% Again loop through cells
for c=1:mesh.TotalCells
    cfaces = mesh.CellFaces{c};
    % Loop through faces in cell
    for ff=1:length(cfaces)
        f = cfaces(ff);
        fid = mesh.FaceID(f);
        % Interior Face
        if fid == 0
            fcells = mesh.FaceCells(f,:);
            if fcells(1) == c
                fnorm = mesh.FaceNormal(f,:)';
                fn1 = dof.FaceCellNodes{f,1};
                fn2 = fliplr(dof.FaceCellNodes{f,2});
                M = fe.FaceMassMatrix{f,1};
            else
                fnorm = - mesh.FaceNormal(f,:)';
                fn1 = dof.FaceCellNodes{f,2};
                fn2 = fliplr(dof.FaceCellNodes{f,1});
                M = fe.FaceMassMatrix{f,2};
            end
        % Boundary Face
        else
            M = fe.FaceMassMatrix{f,1};
            fnorm = mesh.FaceNormal(f,:)';
            fn1 = dof.FaceCellNodes{f,1};
            fn2 = dof.PeriodicFaceDoFs{f}(:,2)';
        end
        % Apply Upwinding Terms by Angle
        for q=1:num_dirs
            fdot = m_quad.AngularDirections(q,:)*fnorm;
            lq1 = fn1;
            % apply mass matrix term based on face normal
            if fdot > 0
                lq2 = lq1;
            else
                lq2 = fn2;
            end
            L{q}(lq1,lq2) = L{q}(lq1,lq2) + fdot*M*PM(lq1,lq1);
        end
    end
end
% Apply angle integration collapse
% ------------------------------------------------------------------------------
for q=1:num_dirs
    T = T + d2m(1,q)*( L{q}\( m2d(1,q)*S ));
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