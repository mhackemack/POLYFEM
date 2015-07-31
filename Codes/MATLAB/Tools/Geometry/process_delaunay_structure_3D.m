%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Tetrahedral Delaunay Structure
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
function obj = process_delaunay_structure_3D(dim,tri)
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
% Get General Triangulation Information
edges = tri.edges();
bfaces = tri.freeBoundary();
bfaces_ord = sort(bfaces,2);
eAttach = tri.edgeAttachments(edges);
tet_ord = tet_face_orderings();
neigh = tri.neighbors();
% Compute Initial Face Information
face_verts = zeros(4*size(cells,1),3);
face_verts_ord = zeros(4*size(cells,1),3);
face_bool = logical(zeros(4*size(cells,1),3));
f = 0;
for c=1:size(cells,1)
    cv = cells(c,:);
    for i=1:4
        f = f + 1;
        fv(f,:) = cv(tet_ord(i,:));
    end
end
fv_ord = sort(fv,2);
[~,fia,fic] = unique(fv_ord,'rows');
% Initialize Initial Data Set
% ---------------------------
obj.TotalCells = size(cells,1);
obj.TotalVertices = size(verts,1);
obj.TotalFaces = length(fia);
obj.TotalBoundaryFaces = size(bfaces,1);
obj.TotalInteriorFaces = obj.TotalFaces - obj.TotalBoundaryFaces;

obj.Vertices = verts;
obj.BodyCenter = mean(obj.Vertices);
obj.MatID = ones(obj.TotalCells,1,'uint32');
obj.FaceID = zeros(obj.TotalFaces,1,'uint32');
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

% Loop through faces and build some information
for f=1:obj.TotalFaces
    obj.FaceVerts{f} = fv(fia(f),:);
    fvo = fv_ord(fia(f),:);
    fverts = verts(obj.FaceVerts{f},:);
    obj.FaceCenter(f,:) = mean(fverts);
    obj.FaceArea(f) = polygonArea3d(fverts);
    for bf = 1:obj.TotalBoundaryFaces
        if fvo == bfaces_ord(bf,:)
            obj.FaceID(f) = 1;
            break;
        end
    end
end

% Loop through cells and build some information
for c=1:obj.TotalCells
    cc = (c-1)*4+1:4*c;
    obj.CellVerts{c} = cells(c,:);
    obj.CellCenter(c,:) = mean(verts(cells(c,:),:));
    obj.CellFaces{c} = fic(cc);
    obj.CellSurfaceArea(c) = sum(obj.FaceArea(obj.CellFaces{c}));
    for f=1:4
        fcc = fic(cc(f));
        if obj.FaceCells(fcc,1) == 0
            obj.FaceCells(fcc,1) = c;
        else
            obj.FaceCells(fcc,2) = c;
        end
        % Form Jacobian and calculate volume
        vv = [obj.Vertices(obj.FaceVerts{fcc},:);obj.CellCenter(c,:)]';
        J = zeros(dim+1);
        J(1:3,:) = vv; J(4,:) = 1;
        obj.CellVolume(c) = obj.CellVolume(c) + abs(det(J))/6;
    end
end

% Calculate Orthogonal Projections
for f=1:obj.TotalFaces
    fcells = obj.FaceCells(f,:);
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

% Calculate Face Normals and Face-Cell Orderings
for f=1:obj.TotalFaces
    fcells = obj.FaceCells(f,:);
    fc = obj.FaceCenter(f,:);
    fv = obj.FaceVerts{f};
    v1 = obj.Vertices(fv(1),:);
    v2 = obj.Vertices(fv(2),:);
    v3 = obj.Vertices(fv(3),:);
    obj.FaceNormal(f,:) = cross(v3-v2, v1-v2);
    if obj.FaceID(f) ~= 0
        cc = obj.CellCenter(fcells(1),:);
        ncell = fc - cc;
        if dot(ncell, obj.FaceNormal(f,:)) < 0
            obj.FaceNormal(f,:) = -1.0*obj.FaceNormal(f,:);
        end
    else
        cc1 = obj.CellCenter(fcells(1),:); fc1 = fc - cc1;
        cc2 = obj.CellCenter(fcells(2),:); fc2 = fc - cc2;
        if dot(fc1, obj.FaceNormal(f,:)) < 0
            obj.FaceNormal(f,:) = -1.0*obj.FaceNormal(f,:);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxiallary Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = tet_face_orderings()
out = [1,2,3;1,3,4;1,4,2;2,4,3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%