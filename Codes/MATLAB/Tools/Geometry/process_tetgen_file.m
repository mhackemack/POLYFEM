function obj = process_tetgen_file(dim,cells,verts,faces,edges)
if nargin < 4
    error('Insufficient Input.')
end
nv = size(verts,2);
if nv == dim
    obj.Vertices = verts;
else
    obj.Vertices = verts(:,1:end-1);
end
% Initialize Initial Data Sets
% ----------------------------
obj.BodyCenter = mean(obj.Vertices);
obj.TotalCells = size(cells,1);
obj.TotalFaces = size(faces,1);
obj.MatID = ones(obj.TotalCells,1,'uint32');
obj.FaceID = uint32(faces(:,end));
obj.CellVerts = cell(obj.TotalCells,1);
obj.CellCenter = zeros(obj.TotalCells,dim);
obj.CellVolume = zeros(obj.TotalCells,1);
obj.CellSurfaceArea = zeros(obj.TotalCells,1);
obj.CellFaces = cell(obj.TotalCells,1);
obj.OrthogonalProjection = zeros(obj.TotalFaces, 2);
obj.FaceVerts = cell(obj.TotalFaces,1);
obj.FaceNormal = zeros(obj.TotalFaces,dim);
obj.FaceCenter = zeros(obj.TotalFaces,dim);
obj.FaceArea = zeros(obj.TotalFaces,1);
obj.FaceCells = zeros(obj.TotalFaces,2,'uint32');
% Loop through cells
for c = 1:obj.TotalCells
    tcell = cells(c,:);
    tvert = obj.Vertices(tcell,:);
    obj.CellVerts{c} = tcell;
    obj.CellCenter(c,:) = mean(tvert);
    [obj.CellVolume(c),obj.CellSurfaceArea(c)] = volume_area_3D(tvert);
end
% Loop through faces/edges
fc = zeros(obj.TotalCells,1);
for f = 1:obj.TotalFaces
    fflag = obj.FaceID(f);
    obj.FaceVerts{f} = faces(f,1:end-1);
    fverts = obj.Vertices(obj.FaceVerts{f},:);
    obj.FaceCenter(f,:) = mean(fverts);
    [~,~,V] = svd(bsxfun(@minus,fverts,mean(fverts,1)));
    obj.FaceNormal(f,:) = V(:,dim)';
    if dot(obj.FaceNormal(f,:),obj.FaceCenter(f,:) - obj.BodyCenter) < 0
        obj.FaceNormal(f,:) = (-1.0)*obj.FaceNormal(f,:);
    end
    obj.FaceArea(f) = polygonArea3d(fverts);
    % Prepare to loop through cells to get face/cell information
    fnodes = (obj.FaceVerts{f});
    nc = 0;
    bool = 0;
    % Loop through cells
    % ------------------
    for c=1:obj.TotalCells
        cnodes = (obj.CellVerts{c});
        if sum(ismember(cnodes,fnodes)) == length(fnodes)
            nc = nc + 1;
            obj.FaceCells(f,nc) = c;
            fc(c) = fc(c) + 1;
            obj.CellFaces{c}(fc(c)) = f;
            if nc == 2
                bool = 1;
            end
        end
        if logical(bool)
            % If interior cell - reorient cells so that FaceCells(1)
            % is the "positive" cell in relation to the face normal.
            if fflag == 0
                fnorm = obj.FaceNormal(f,:);
                fcenter = obj.FaceCenter(f,:);
                ccenter = obj.CellCenter(obj.FaceCells(f,:),:);
                e = dot(fnorm, ccenter(1,:) - fcenter);
                if e < 0
                    obj.FaceCells(f,:) = fliplr(obj.FaceCells(f,:));
                end
            end
            break
        end
    end
    % Orthogonal Projections - used for MIP
    fcells = obj.FaceCells(f,:);
    if fflag == 0
        for c=1:2
            nv = length(obj.CellVerts{fcells(c)});
            obj.OrthogonalProjection(f,c) = get_orthogonal_length(dim, ...
                obj.CellVolume(fcells(c)), nv, obj.FaceArea(f), obj.CellSurfaceArea(fcells(c)));
        end
    else
        nv = length(obj.CellVerts{fcells(1)});
        obj.OrthogonalProjection(f,1) = get_orthogonal_length(dim, ...
                obj.CellVolume(fcells(1)), nv, obj.FaceArea(f), obj.CellSurfaceArea(fcells(1)));
    end
end

