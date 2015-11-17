%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Polymesher Mesh Input
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB function to read the output of the PolyMesher suite
%                   and convert to our mesh structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj  = process_polymesher_structure(N, E)
obj.Dimension = 2;
% Get Number of Vertices and Elements
obj.TotalVertices = size(N,1);
obj.TotalCells = length(E);
vertbools = false(obj.TotalVertices,1);
obj.Vertices = N;
xmin = min(N(:,1)); xmax = max(N(:,1));
ymin = min(N(:,2)); ymax = max(N(:,2));
% Allocate Cell Arrays
obj.MatID = ones(obj.TotalCells,1,'uint32');
obj.CellVerts = cell(obj.TotalCells,1);
obj.CellCenter = zeros(obj.TotalCells,obj.Dimension);
obj.CellVolume = zeros(obj.TotalCells,1);
obj.CellFaces = cell(obj.TotalCells,1);
obj.CellFaceVerts = cell(obj.TotalCells,1);
obj.CellSurfaceArea = zeros(obj.TotalCells,1);
obj.CellFaces = cell(obj.TotalCells,1);
% Loop through cells and get info - also get temp face verts array
t_face_verts = [];
t_face_verts_ord = [];
t_face_cells = [];
for c=1:obj.TotalCells
    vertbools(E{c}) = true;
    obj.CellVerts{c} = E{c};
    obj.CellCenter(c,:) = mean( N(E{c},:) );
    obj.CellVolume(c) = polygonArea(obj.Vertices(obj.CellVerts{c},:));
    % Loop through cell verts to build temp face array
    nv = length(E{c});
    for i=1:length(E{c})
        ii = [i,mod(i,nv)+1];
        fverts = E{c}(ii);
        t_face_verts = [t_face_verts; fverts];
        t_face_cells = [t_face_cells; c];
        if fverts(1) < fverts(2)
            t_face_verts_ord = [t_face_verts_ord; [fverts(1),fverts(2)]];
        else
            t_face_verts_ord = [t_face_verts_ord; [fverts(2),fverts(1)]];
        end
    end
end
% Get unique ordering of face information
[C,IA,IC] = unique(t_face_verts_ord, 'stable', 'rows');
obj.TotalFaces = size(C,1);
% Generate Face Arrays
obj.FaceVerts = cell(obj.TotalFaces,1);
obj.FaceNormal = zeros(obj.TotalFaces,obj.Dimension);
obj.FaceCenter = zeros(obj.TotalFaces,obj.Dimension);
obj.FaceArea = zeros(obj.TotalFaces,1);
obj.FaceCells = zeros(obj.TotalFaces,2);
obj.FaceID = zeros(obj.TotalFaces,1,'uint32');
obj.OrthogonalProjection = zeros(obj.TotalFaces, 2);
% Loop through faces and assign face info
for f=1:obj.TotalFaces
    iia = t_face_verts(IA(f),:);
    obj.FaceVerts{f} = iia;
    obj.FaceCenter(f,:) = mean( N(iia,:) ); fc = obj.FaceCenter(f,:);
    fvv = obj.Vertices(obj.FaceVerts{f},:); dx = diff(fvv);
    obj.FaceArea(f) = norm(dx);
    obj.FaceNormal(f,:) = [dx(2),-dx(1)];
    obj.FaceNormal(f,:) = obj.FaceNormal(f,:) / norm(obj.FaceNormal(f,:));
    if abs(xmin - fc(1)) < 1e-14 || abs(xmax - fc(1)) < 1e-14 || ...
       abs(ymin - fc(2)) < 1e-14 || abs(ymax - fc(2)) < 1e-14
        obj.FaceID(f) = 1;
    end
end
% Loop through cells again to assign faces 
cf_count = 0;
for c=1:obj.TotalCells
    % Loop through cell faces
    nv = length(obj.CellVerts{c});
    cft = []; sa = 0;
    for f=1:length(obj.CellVerts{c})
        cf_count = cf_count + 1;
        ff = IC(cf_count);
        cft = [cft,ff];
        sa = sa + obj.FaceArea(ff);
        if obj.FaceCells(ff,1) == 0
            obj.FaceCells(ff,1) = c;
        else
            obj.FaceCells(ff,2) = c;
        end
    end
    obj.CellFaces{c} = cft;
    obj.CellSurfaceArea(c) = sa;
    % Loop through faces again for orthogonal projection
    for f=1:length(obj.CellFaces{c})
        ff = obj.CellFaces{c}(f);
        h = get_orthogonal_length(2,obj.CellVolume(c),nv,obj.FaceArea(ff),sa);
        if obj.FaceCells(ff,1) == c
            obj.OrthogonalProjection(ff,1) = h;
        elseif obj.FaceCells(ff,2) == c
            obj.OrthogonalProjection(ff,2) = h;
        end
    end
end
% Loop through faces 1 last time to check normal directions
for f=1:obj.TotalFaces
    c1 = obj.FaceCells(f,1);
    fc = obj.FaceCenter(f,:);
    cc = obj.CellCenter(c1,:);
    dx = fc - cc;
    if dot(dx, obj.FaceNormal(f,:)) < 0
        obj.FaceNormal(f,:) = -1.0*obj.FaceNormal(f,:);
    end
end