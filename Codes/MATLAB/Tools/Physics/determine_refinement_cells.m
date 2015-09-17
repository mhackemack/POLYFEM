%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate Refinement Cells
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to determine which cells should be refined
%                   based on some method.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function determine_refinement_cells(data, mesh, DoF, FE, x)
global glob
% Process Inputs
if ~iscell(x)
    y = cell(1,1);
    y{1} = x;
else
    y = x;
end
ny = length(y);
clear x
if strcmp(data.Neutronics.transportMethod,'Diffusion')
    diff_bool = true;
    D = data.Neutronics.Diffusion.DiffXS;
elseif strcmp(data.Neutronics.transportMethod,'Transport')
    D = ones(data.problem.NumberMaterials,ny);
    diff_bool = false;
end
err = zeros(mesh.TotalCells,1);
% Determine Mesh Refinement Based on FEM Type
%   CFEM = Gradient Jumps on Edges
%   DFEM = Solution Jumps on Edges
% -------------------------------------------
if DoF.FEMType == 1
    face_diff = zeros(mesh.TotalFaces,1);
    % Loop through faces in mesh
    % --------------------------
    for f=1:mesh.TotalFaces
        flag = mesh.FaceID(f);
        fnorm = mesh.FaceNormal(f,:);
        c1 = mesh.FaceCells(f,1); m1 = mesh.MatID(c1);
        if flag == 0
            c2 = mesh.FaceCells(f,2);       m2 = mesh.MatID(c2);
            q1 = FE.FaceQuadNodes{f,1};     q2 = FE.FaceQuadNodes{f,2};
            w1 = FE.FaceQuadWeights{f,1};   w2 = FE.FaceQuadWeights{f,2};
            g1 = FE.FaceBasisGrads{f,1};    g2 = FE.FaceBasisGrads{f,2};
            % Determine across face quadrature order - 2D is simply the opposite
            % orientation. 3D is more complicated if I ever do it.
            if mesh.Dimension == 1
                q_ind = 1;
            elseif mesh.Dimension == 2
                q_ind = size(q2,1):-1:1;
            else
                
            end
            % Loop through nodes and compute estimator
            for q=1:size(q1,1)
                cn1 = DoF.ConnectivityArray{c1}; cn2 = DoF.ConnectivityArray{c2};
                w = w1(q);
                tg1 = g1(:,:,q); tg2 = g2(:,:,q_ind(q));
                for i=1:ny
                    gy1 = y{i,1}(cn1)'*tg1; gy2 = y{i,1}(cn2)'*tg2;
                    dg = fnorm*(D(m1,i)*gy1-D(m2,i)*gy2)';
                    face_diff(f) = face_diff(f) + w*dg^2;
                end
            end
        elseif diff_bool
            if data.Neutronics.Diffusion.BCFlags(flag) ~= glob.Neumann, continue; end
            q = FE.FaceQuadNodes{f,1};
            w = FE.FaceQuadWeights{f,1};
            g = FE.FaceBasisGrads{f,1};
            cn = DoF.ConnectivityArray{c1};
            for i=1:ny
                
            end
        end
    end
    % Loop through cells in mesh
    for c=1:mesh.TotalCells
        faces = mesh.CellFaces{c};
        % Loop through faces in cell
        % --------------------------
        for ff=1:length(faces)
            f = faces(ff);
            if c == mesh.FaceCells(f,1)
                err(c) = err(c) + sqrt(mesh.OrthogonalProjection(f,1)/24*face_diff(f));
            elseif c == mesh.FaceCells(f,2)
                err(c) = err(c) + sqrt(mesh.OrthogonalProjection(f,2)/24*face_diff(f));
            end
        end
    end
elseif DoF.FEMType == 2
    % Loop through cells in mesh - estimate error for each cell
    % ---------------------------------------------------------
    for c=1:mesh.TotalCells
        faces = mesh.CellFaces{c};
        % Loop through faces in cell
        % --------------------------
        for ff=1:length(faces)
            f = faces(ff);
            flag = mesh.FaceID(f);
            % Only perform measurement over interior faces
            if flag == 0
                ndp = DoF.FaceNodePartners{f};
                if c == mesh.FaceCells(f,1)
                    M = FE.FaceMassMatrix{f,1};
                elseif c == mesh.FaceCells(f,2)
                    M = FE.FaceMassMatrix{f,2};
                end
                for iy=1:ny
                    err(c) = err(c) + sum(M*abs(y{iy,1}(ndp(:,1)) - y{iy,1}(ndp(:,2))));
                end
            end
        end
        err(c) = abs(err(c)) / abs(mesh.CellSurfaceArea(c));
    end
end
% Normal error estimates and set refinement flags
err = err / max(err);
% Switch based on cell refinement type
if data.problem.refinementType == 0
    % Loop through cells in mesh
    for c=1:mesh.TotalCells
        if err(c) >= data.problem.refinementTolerance
            mesh.set_refinement_flag(c);
        end
    end
elseif data.problem.refinementType == 1
    % Sort cell errors
    [~,ind] = sort(err,1,'descend');
    cell_list = (1:mesh.TotalCells)';
    cell_list = cell_list(ind);
    num = round((data.problem.refinementTolerance)*mesh.TotalCells);
    for c=1:num
        mesh.set_refinement_flag(cell_list(c));
    end
end