function obj = process_triangle_file(dim,cells,verts,faces)
if nargin < 3
    error('Insufficient Input.')
end
obj.Vertices = verts(:,1:end-1);
obj.BodyCenter = mean(obj.Vertices);
obj.TotalCells = size(cells,1);
obj.TotalFaces = size(faces,1);
obj.MatID = ones(obj.TotalCells,1,'uint32');
obj.FaceID = uint32(faces(:,end));
obj.CellVerts = cell(obj.TotalCells,1);
obj.CellCenter = zeros(obj.TotalCells,2);
obj.CellVolume = zeros(obj.TotalCells,1);
obj.CellFaceVerts = cell(obj.TotalCells,1);
obj.CellSurfaceArea = zeros(obj.TotalCells,1);
obj.CellFaces = cell(obj.TotalCells,1);
obj.OrthogonalProjection = zeros(obj.TotalFaces, 2);
obj.FaceVerts = cell(obj.TotalFaces,1);
obj.FaceNormal = zeros(obj.TotalFaces,2);
obj.FaceCenter = zeros(obj.TotalFaces,2);
obj.FaceArea = zeros(obj.TotalFaces,1);
obj.FaceCells = zeros(obj.TotalFaces,2,'uint32');
% Loop through cells
for c = 1:obj.TotalCells
    tcell = cells(c,:);
    tvert = obj.Vertices(tcell,:);
    [~,ix] = sort(atan2(tvert(:,2)-mean(tvert(:,2)),tvert(:,1)-mean(tvert(:,1))));
    obj.CellVerts{c} = tcell(ix');
    obj.CellCenter(c,:) = mean(tvert);
    obj.CellVolume(c) = triangle_area(tvert);
    obj.CellFaceVerts{c}{1} = obj.CellVerts{c}([1,2]);
    obj.CellFaceVerts{c}{2} = obj.CellVerts{c}([2,3]);
    obj.CellFaceVerts{c}{3} = obj.CellVerts{c}([3,1]);
end
% Loop through faces/edges
for f = 1:obj.TotalFaces
    obj.FaceVerts{f} = faces(f,1:end-1);
    obj.FaceCenter(f,:) = mean(obj.Vertices(obj.FaceVerts{f},:));
    dd = obj.Vertices(obj.FaceVerts{f}(2),:)-obj.Vertices(obj.FaceVerts{f}(1),:);
    obj.FaceNormal(f,:) = [dd(2),-dd(1)];
    obj.FaceNormal(f,:) = obj.FaceNormal(f,:) / norm(obj.FaceNormal(f,:));
    if dot(obj.FaceNormal(f,:),obj.FaceCenter(f,:) - obj.BodyCenter) < 0
        obj.FaceNormal(f,:) = (-1.0)*obj.FaceNormal(f,:);
    end
    obj.FaceArea(f) = norm(dd);
end
% Loop through faces again
fc = zeros(obj.TotalCells,1);
for f = 1:obj.TotalFaces
    fflag = obj.FaceID(f);
    nodes = sort(obj.FaceVerts{f});
    nc = 0;
    bool = 0;
    for c=1:obj.TotalCells
        for v=1:3
            if v==3
                vv = obj.CellVerts{c}([3,1]);
            else
                vv = obj.CellVerts{c}([v,v+1]);
            end
            if vv(1) > vv(2)
                vtemp = vv(1);
                vv(1) = vv(2);
                vv(2) = vtemp;
            end
            if vv(1) == nodes(1) && vv(2) == nodes(2)
                nc = nc + 1;
                obj.FaceCells(f,nc) = c;
                fc(c) = fc(c) + 1;
                obj.CellFaces{c}(fc(c)) = f;
                if nc == 2
                    bool = 1;
                end
            end
        end
        if logical(bool)
            % If interior cell - reorient cells so that FaceCells(1)
            % is the "positive" cell in relation to the face normal.
            if fflag == 0
                fnorm = obj.FaceNormal(f,:);
                fcenter = obj.FaceCenter(f,:);
                ccenter = obj.CellCenter(obj.FaceCells(f,:),:);
                e = dot(fnorm, fcenter - ccenter(1,:));
                if e < 0
                    obj.FaceCells(f,:) = fliplr(obj.FaceCells(f,:));
                end
            end
            break
        end
    end
end
% Loop through cells
for c=1:obj.TotalCells
    cfaces = obj.CellFaces{c};
    for f=1:length(cfaces)
        obj.CellSurfaceArea(c) = obj.CellSurfaceArea(c) + obj.FaceArea(cfaces(f));
    end
end
% Loop through faces/edges
for f = 1:obj.TotalFaces
    fflag = obj.FaceID(f);
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

function out = triangle_area(verts)
% This is Heron's formula
a = norm(verts(2,:) - verts(1,:));
b = norm(verts(3,:) - verts(2,:));
c = norm(verts(1,:) - verts(3,:));
s = (a+b+c)/2;
out = sqrt(s*(s-c)*(s-b)*(s-a));