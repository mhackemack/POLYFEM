%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execution Functor - LHS DFEM Transport (upwind)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = exec_func_LHS_dfem_transport_upwind(data, mesh, DoF, FE, angs, groups)
% Process Input Space
% -------------------
global glob
dim = mesh.Dimension;
ndat = data.Transport;
ndof = DoF.TotalDoFs;
ng = length(groups);
na = length(angs);
ntot = ndof * ng * na;
angdirs = ndat.AngularDirections;
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
% Allocate Memory Space
% ---------------------
if ntot > glob.maxMatrix
    L = get_sparse_matrix(data, mesh, DoF, FE, angs, groups);
    return
end
L = zeros(ntot);
% Loop through cells
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    cnodes = DoF.ConnectivityArray{c};
    M = FE.CellMassMatrix{c};
    G = FE.CellGradientMatrix{c};
    % Loop through angles
    for q=1:na
        tq = angs(q);
        GG = cell_dot(dim, angdirs(tq, :), G);
        % Loop through energy groups
        for g=1:ng
            cnqg = cnodes + g_offset(g) + q_offset(q);
            L(cnqg,cnqg) = L(cnqg,cnqg) + ndat.TotalXS(cmat,groups(g))*M + GG;
        end
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
        for qq=1:na
            sxs = ndat.ScatteringXS(cmat,1,1,1)*ndat.discrete_to_moment(qq);
            ccc = cnodes + q_offset(qq);
            L(cnqg,ccc) = L(cnqg,ccc) - sxs*M*ndat.moment_to_discrete(tq);
        end
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
    end
end

% Loop through faces
for f=1:mesh.TotalFaces
    fflag = mesh.FaceID(f);
    fnorm = mesh.FaceNormal(f,:);
    % Perform Interior Face operations
    if fflag == 0
        % Loop through angles
        for q=1:na
            tq = angs(q);
            fdot = fnorm * angdirs(tq,:)';
            if fdot < 0
                fdot = -1*fdot;
                fnodes1 = DoF.FaceCellNodes{f,1};
                fnodes2 = DoF.FaceCellNodes{f,2};
                M = FE.FaceMassMatrix{f,1};
                MM = FE.FaceConformingMassMatrix{f,1};
            else
                fnodes1 = DoF.FaceCellNodes{f,2};
                fnodes2 = DoF.FaceCellNodes{f,1};
                M = FE.FaceMassMatrix{f,2};
                MM = FE.FaceConformingMassMatrix{f,2};
            end
            for g=1:ng
                fnqg1 = fnodes1 + g_offset(g) + q_offset(q);
                fnqg2 = fnodes2 + g_offset(g) + q_offset(q);
                L(fnqg1,fnqg1) = L(fnqg1,fnqg1) + fdot*M;
                L(fnqg1,fnqg2) = L(fnqg1,fnqg2) - fdot*MM;
            end
        end
    % Perform Boundary Face operations
    else
        fnorm = mesh.FaceNormal(f,:);
        fnodes = DoF.FaceCellNodes{f,1};
        M = FE.FaceMassMatrix{f,1};
        % Loop through angles
        for q=1:na
            tq = angs(q);
            fdot = fnorm * angdirs(tq,:)';
            if fdot < 0 && abs(fdot) > glob.small
                % Loop through energy groups
                for g=1:ng
                    fnqg = fnodes + g_offset(g) + q_offset(q);
                    L(fnqg,fnqg) = L(fnqg,fnqg) - M*fdot;
                end
            end
        end
    end
end
if ~issparse(L), L = sparse(L); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = get_sparse_matrix(data, mesh, DoF, FE, angs, groups)
% Process Input Space
% -------------------
global glob
dim = mesh.Dimension;
ndat = data.Transport;
ndof = DoF.TotalDoFs;
ng = length(groups);
na = length(angs);
ntot = ndof * ng * na;
angdirs = ndat.AngularDirections;
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
% Allocate Memory Space
% ---------------------
I = [];
J = [];
TMAT = [];
% Loop through cells
for c=1:mesh.TotalCells
    cmat = mesh.MatID(c);
    cnodes = DoF.ConnectivityArray{c}; ncnodes = length(cnodes);
    onesnodes = ones(ncnodes,1);
    M = FE.CellMassMatrix{c};
    G = FE.CellGradientMatrix{c};
    % Loop through angles
    for q=1:na
        tq = angs(q);
        GG = cell_dot(dim, angdirs(tq, :), G);
        % Loop through energy groups
        for g=1:ng
            % Get Indexing Information
            % ------------------------
            cnqg = cnodes + g_offset(g) + q_offset(q);
            rows = (onesnodes*cnqg)';
            cols = onesnodes*cnqg;
            % Apply Sparse Indexing
            % ---------------------
            tmat = ndat.TotalXS(cmat,g)*M + GG;
            I = [I;rows(:)];
            J = [J;cols(:)];
            TMAT = [TMAT;tmat(:)];
        end
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
%         for qq=1:na
%             tq = angs(qq);
%             ccc = cnodes + q_offset(qq);
%             cols2 = onesnodes*ccc;
%             I = [I;rows(:)];
%             J = [J;cols2(:)];
%             tmat = -ndat.ScatteringXS(cmat,1,1,1)*ndat.discrete_to_moment(tq)*M*ndat.moment_to_discrete(tq);
%             TMAT = [TMAT;tmat(:)];
%         end
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
        % UNCOMMENT THIS FOR MONOCHROMATIC SCATTERING!!!
    end
end
% Loop through faces
for f=1:mesh.TotalFaces
    fflag = mesh.FaceID(f);
    fnorm = mesh.FaceNormal(f,:);
    % Perform Interior Face operations
    if fflag == 0
        % Loop through angles
        for q=1:na
            tq = angs(q);
            fdot = fnorm * angdirs(tq,:)';
            if fdot < 0
                fdot = -1*fdot;
                fnodes1 = DoF.FaceCellNodes{f,1};
                fnodes2 = DoF.FaceCellNodes{f,2};
                M = FE.FaceMassMatrix{f,1};
                MM = FE.FaceConformingMassMatrix{f,1};
            else
                fnodes1 = DoF.FaceCellNodes{f,2};
                fnodes2 = DoF.FaceCellNodes{f,1};
                M = FE.FaceMassMatrix{f,2};
                MM = FE.FaceConformingMassMatrix{f,2};
            end
            onesnodes = ones(length(fnodes1),1);
            for g=1:ng
                fnqg1 = fnodes1 + g_offset(g) + q_offset(q);
                fnqg2 = fnodes2 + g_offset(g) + q_offset(q);
                % Apply within cell term
                % ----------------------
                cols1 = onesnodes*fnqg1; rows1 = (onesnodes*fnqg1)'; tmat = fdot*M;
                I = [I;rows1(:)]; J = [J;cols1(:)]; TMAT = [TMAT;tmat(:)];
                % Apply cross-cell term
                % ---------------------
                cols2 = onesnodes*fnqg2; tmat = -1.0*fdot*MM;
                I = [I;rows1(:)]; J = [J;cols2(:)]; TMAT = [TMAT;tmat(:)];
            end
        end
    % Perform Boundary Face operations
    else
        fnorm = mesh.FaceNormal(f,:);
        fnodes = DoF.FaceCellNodes{f,1};
        onesnodes = ones(length(fnodes),1);
        M = FE.FaceMassMatrix{f,1};
        % Loop through angles
        for q=1:na
            tq = angs(q);
            fdot = fnorm * angdirs(tq,:)';
            if fdot < 0 && abs(fdot) > glob.small
                % Loop through energy groups
                for g=1:ng
                    fnqg = fnodes + g_offset(g) + q_offset(q);
                    cols = onesnodes*fnqg; rows = (onesnodes*fnqg)'; tmat = -fdot*M;
                    I = [I;rows(:)]; J = [J;cols(:)]; TMAT = [TMAT;tmat(:)];
                end
            end
        end
    end
end
L = sparse(I,J,TMAT,ntot,ntot);
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