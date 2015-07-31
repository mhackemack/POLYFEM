%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Triangle Delaunay Structure
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate all data structures necessary
%                   to fully describe a geometric domain to be used for
%                   finite element (FEM) calculations.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = process_delaunay_structure_2D(dim,tri)
if isa(tri, 'delaunayTriangulation')
    verts = tri.Points;
    cells = tri.ConnectivityList;
elseif isa(tri, 'DelaunayTri')
    verts = tri.X;
    cells = tri.Triangulation;
elseif isa(tri, 'triangulation')
    verts = tri.Points;
    cells = tri.ConnectivityList;
else
    error('Could not determine delaunay structure type.')
end
% Get Face Information
faces = tri.edges();
bfaces = tri.freeBoundary();
eAttach = tri.edgeAttachments(faces);
% Initialize Initial Data Set
% ---------------------------
obj.TotalCells = size(cells,1);
obj.TotalFaces = size(faces,1);
obj.TotalVertices = size(verts,1);

obj.Vertices = verts;
obj.BodyCenter = mean(obj.Vertices);
obj.MatID = ones(obj.TotalCells,1,'uint32');
obj.FaceID = zeros(obj.TotalFaces,1,'uint32');
obj.CellVerts = cell(obj.TotalCells,1);
obj.CellCenter = zeros(obj.TotalCells,dim);
obj.CellVolume = zeros(obj.TotalCells,1);
obj.CellFaceVerts = cell(obj.TotalCells,1);
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
    [obj.CellVolume(c),obj.CellSurfaceArea(c)] = polyareaperimeter(tvert);
    obj.CellFaceVerts{c}{1} = obj.CellVerts{c}([1,2]);
    obj.CellFaceVerts{c}{2} = obj.CellVerts{c}([2,3]);
    obj.CellFaceVerts{c}{3} = obj.CellVerts{c}([3,1]);
end
obj.CellVolume = abs(obj.CellVolume);
% Loop through faces/edges
for f = 1:obj.TotalFaces
    if dim == 2
        ffaces = faces(f,:);
        if ffaces(2) > ffaces(1)
            t = ffaces(2);
            ffaces(2) = ffaces(1);
            ffaces(1) = t;
        end
        obj.FaceVerts{f} = faces(f,:);
        fverts = obj.Vertices(obj.FaceVerts{f},:);
        for i=1:size(bfaces,1)
            bbfaces = bfaces(i,:);
            if bbfaces(2) > bbfaces(1)
                t = bbfaces(2);
                bbfaces(2) = bbfaces(1);
                bbfaces(1) = t;
            end
            if ffaces == bbfaces
                obj.FaceID(f) = 1;
                bfaces(i,:) = [];
                break
            end
        end
    else
        ffaces = sort(faces(f,:));
        obj.FaceVerts{f} = faces(f,:);
        fverts = obj.Vertices(obj.FaceVerts{f},:);
        for i=1:length(bfaces)
            if ffaces == sort(bfaces(i,:))
                obj.FaceID(f) = 1;
                bfaces(i,:) = [];
                break
            end
        end
    end
    obj.FaceCenter(f,:) = mean(fverts);
    [~,~,V] = svd(bsxfun(@minus,fverts,mean(fverts,1)));
    obj.FaceNormal(f,:) = V(:,dim)';
    if dim == 2
        obj.FaceArea(f) = norm(diff(fverts));
    elseif dim == 3
        obj.FaceArea(f) = polygonArea3d(fverts);
    end
    obj.FaceCells(f,:) = eAttach{f};
    for c=1:length(eAttach{f})
        obj.CellFaces{eAttach{f}(c)} = [obj.CellFaces{eAttach{f}(c)},f];
    end
    % Check normal ordering
    fcells = obj.FaceCells(f,:);
    fnorm = obj.FaceNormal(f,:);
    tnorm = (obj.FaceCenter(f,:) - obj.CellCenter(fcells(1),:)); tnorm = tnorm / norm(tnorm);
    if dot(fnorm, tnorm) < 0
        obj.FaceNormal(f,:) = -1.0*obj.FaceNormal(f,:);
    end
    % Orthogonal Projections - used for MIP
    if obj.FaceID(f) == 0
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
% perform cleanup
for f = 1:obj.TotalFaces
    tcell = obj.FaceCells(f,1);
    fverts = obj.FaceVerts{f};
    cverts = obj.CellVerts{tcell};
    for i=1:length(cverts)
        if cverts(i) == fverts(1)
            if i==length(cverts)
                if cverts(1) ~= fverts(2)
                    obj.FaceVerts{f} = fliplr(obj.FaceVerts{f});
                end
            else
                if cverts(i+1) ~= fverts(2)
                    obj.FaceVerts{f} = fliplr(obj.FaceVerts{f});
                end
            end
        end
    end
end

