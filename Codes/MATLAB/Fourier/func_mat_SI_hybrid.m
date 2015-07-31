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
function [L, S] = func_mat_SI_hybrid(lam, input)
% Copy Input Space
% ----------------
data = input.data;
mesh = input.mesh;
dof = input.dof;
fe = input.fe;
m_quad = input.Quadrature;
m2d = m_quad.moment_to_discrete;
d2m = m_quad.discrete_to_moment;
angnorm = m_quad.AngQuadNorm;
g0 = data.FluxStabilization;
d0 = data.CurrentStabilization;
% Retrieve Preliminary Data
% -------------------------
diamD = mesh.Diameter;
dim = mesh.Dimension;
ndofs = dof.TotalDoFs;
node_locs = dof.NodeLocations;
if dim == size(lam,2); lam=lam'; end
num_dirs = m_quad.NumberAngularDirections;
q_offset = (1:num_dirs)*ndofs - ndofs;
PV = exp(1i*node_locs*lam);
PM = diag(PV);
% Allocate Matrix Arrays
% ----------------------
ntot = ndofs*num_dirs;
S = zeros(ntot);
L = zeros(ntot);
% Loop through Cells and Build Volumetric Portions of Matrices
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cn  = dof.ConnectivityArray{c};
    mat = mesh.MatID(c);
    M   = fe.CellMassMatrix{c};
    G   = fe.CellGradientMatrix{c};
    txs = data.TotalXS(mat);
    sxs = data.ScatteringXS(mat);
    % Loop through quadrature
    for q=1:num_dirs
        cnq = cn + q_offset(q);
        tm2d = m2d(1,q);
        GG = cell_dot(dim, m_quad.AngularDirections(q,:), G);
%         L(cnq,cnq) = L(cnq,cnq) + (txs*M - GG')*PM(cn,cn);
        L(cnq,cnq) = L(cnq,cnq) + (txs*M + GG)*PM(cn,cn);
        for qq=1:num_dirs
            cnqq = cn + q_offset(qq);
            S(cnq,cnqq) = S(cnq,cnqq) + tm2d*sxs*d2m(1,qq)*M*PM(cn,cn);
        end
    end
end
% Build Face Contributions to Transport Matrix
% ------------------------------------------------------------------------------
for f=1:mesh.TotalFaces
    fid = mesh.FaceID(f);
    fnorm = mesh.FaceNormal(f,:)';
    fcells = mesh.FaceCells(f,:);
    % Interior Face
    if fid == 0
        mID    = mesh.MatID(fcells);
        fn1    = dof.FaceCellNodes{f,1};
        fn2    = fliplr(dof.FaceCellNodes{f,2});
        h      = mesh.OrthogonalProjection(f,:);
        M      = fe.FaceMassMatrix{f,1};
    % Boundary Face
    else
        op_f = mesh.PeriodicOppositeFaces(f);
        op_c = mesh.FaceCells(op_f,1);
        mID  = mesh.MatID([fcells(1), op_c]);
        fn1  = dof.FaceCellNodes{f,1};
        fn2  = dof.PeriodicFaceDoFs{f}(:,2)';
        h    = mesh.OrthogonalProjection([f,op_f],1);
        M    = fe.FaceMassMatrix{f,1};
    end
    if data.StabilizationType == 0
        gm = 1; gp = 1;
        dm = 0; dp = 0;
    elseif data.StabilizationType == 1
        sigs1 = data.ScatteringXS(mID(1));
        sigs2 = data.ScatteringXS(mID(2));
        gm = g0/max(g0,sigs1*diamD);
        gp = g0/max(g0,sigs2*diamD);
        dm = d0*(1-gm)/gm; dp = d0*(1-gp)/gp;
    elseif data.StabilizationType == 2
        sigs1 = data.ScatteringXS(mID(1));
        sigs2 = data.ScatteringXS(mID(2));
        gm = g0/max(g0,sigs1*h(1)); 
        gp = g0/max(g0,sigs2*h(2));
        dm = 0; dp = 0;
    elseif data.StabilizationType == 3
        sigs1 = data.ScatteringXS(mID(1));
        sigs2 = data.ScatteringXS(mID(2));
        gm = g0/max(g0,sigs1*h(1)); 
        gp = g0/max(g0,sigs2*h(2));
        dm = d0*(1-gm)/gm; dp = d0*(1-gp)/gp;
    end
    % Apply Upwinding Terms by Angle
    for q=1:num_dirs
        fnq1 = fn1 + q_offset(q);
        fnq2 = fn2 + q_offset(q);
        fdot = m_quad.AngularDirections(q,:)*fnorm; afdot = abs(fdot);
        PPM = PM(fn1,fn1);
        
        
        
%         % [[u]] terms
%         L(fnq1,fnq1) = L(fnq1,fnq1) - gm*afdot/2*M*PPM; % (-,-)
%         L(fnq1,fnq2) = L(fnq1,fnq2) + gp*afdot/2*M*PPM; % (-,+)
%         % {{u}} terms
%         L(fnq1,fnq1) = L(fnq1,fnq1) + fdot/2*M*PPM; % (-,-)
%         L(fnq1,fnq2) = L(fnq1,fnq2) + fdot/2*M*PPM; % (-,+)
%         % {{Jun}} terms
%         for qq=1:num_dirs
%             fnqq1 = fn1 + q_offset(qq);
%             fnqq2 = fn2 + q_offset(qq);
%             ffdot = m_quad.AngularDirections(qq,:)*fnorm;
%             wt = m_quad.AngularWeights(qq);
%             L(fnq1,fnqq1) = L(fnq1,fnqq1) + dm*fdot*ffdot*wt*M/angnorm*PPM; % (-,-)
%             L(fnq1,fnqq2) = L(fnq1,fnqq2) - dp*fdot*ffdot*wt*M/angnorm*PPM; % (-,+)
%         end
        
        
        
        % ( [[u]] , [[b]] ) terms
        L(fnq1,fnq1) = L(fnq1,fnq1) + gm*afdot/2*M*PPM; % (-,-)
        L(fnq1,fnq2) = L(fnq1,fnq2) - gp*afdot/2*M*PPM; % (-,+)
        % ( {{un}} , {{b}} ) terms
        L(fnq1,fnq1) = L(fnq1,fnq1) - fdot/2*M*PPM; % (-,-)
        L(fnq1,fnq2) = L(fnq1,fnq2) + fdot/2*M*PPM; % (-,+)
        % ( {{Jun}} , {{Jbn}} ) terms
        for qq=1:num_dirs
            fnqq1 = fn1 + q_offset(qq);
            fnqq2 = fn2 + q_offset(qq);
            ffdot = m_quad.AngularDirections(qq,:)*fnorm;
            wt = m_quad.AngularWeights(qq);
            L(fnq1,fnqq1) = L(fnq1,fnqq1) + dm*fdot*ffdot*wt*M/angnorm*PPM; % (-,-)
            L(fnq1,fnqq2) = L(fnq1,fnqq2) - dp*fdot*ffdot*wt*M/angnorm*PPM; % (-,+)
        end
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