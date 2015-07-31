%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process 2D Discontinuous Galerkin Mesh Input
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the Maximum-Entropy basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = process_2D_DG_structure(s_name)
obj.Dimension = 2;
% Get File Text
fid=fopen(s_name);
ttext = textscan(fid,'%n','commentstyle','#'); 
ftext = ttext{1}; clear ttext;
% Get Top of File Info
xmin = 0; ymin = 0;
xmax = ftext(1);
ymax = ftext(2);
obj.TotalCells = ftext(3);
% Initialize Cell Arrays
cell_dg_verts = cell(obj.TotalCells, 1);
n_vertices = zeros(obj.TotalCells, 1);
obj.MatID = ones(obj.TotalCells,1,'uint32');
obj.CellVerts = cell(obj.TotalCells,1);
obj.CellCenter = zeros(obj.TotalCells,obj.Dimension);
obj.CellVolume = zeros(obj.TotalCells,1);
obj.CellFaceVerts = cell(obj.TotalCells,1);
obj.CellSurfaceArea = zeros(obj.TotalCells,1);
obj.CellFaces = cell(obj.TotalCells,1);
% Read DG Connectivity
ind = 4;
for c=1:obj.TotalCells
    n_vertices(c)=ftext(ind);
    i1=ind+1;
    i2=ind+n_vertices(c);
    cell_dg_verts{c} = zeros(n_vertices(c),1);
    cell_dg_verts{c}(:)=ftext(i1:i2);
    ind=i2+1;
    obj.MatID(c) = ftext(ind);
    ind = ind + 2;
end
% Read DG Vertices
n_vert = ftext(ind);
vert = zeros(n_vert, 2);
ind = ind + 1;
for i=1:n_vert
    vert(i,1)=ftext(ind);
    ind = ind + 1;
    vert(i,2)=ftext(ind);
    ind = ind + 1;
end
% Read Mesh Vertices
obj.TotalVertices = (length(ftext)-ind)/2;
% obj.TotalVertices = ftext(ind);
obj.Vertices = zeros(obj.TotalVertices, obj.Dimension);
obj.BodyCenter = [(xmin+xmax)/2, (ymin+ymax)/2];
ind = ind + 1;
for i=1:obj.TotalVertices
    obj.Vertices(i,1)=ftext(ind);
    ind = ind + 1;
    obj.Vertices(i,2)=ftext(ind);
    ind = ind + 1;
    if abs(obj.Vertices(i,1) - xmin) < 1e-14
        obj.Vertices(i,1) = xmin;
    end
    if abs(obj.Vertices(i,1) - xmax) < 1e-14
        obj.Vertices(i,1) = xmax;
    end
    if abs(obj.Vertices(i,2) - ymin) < 1e-14
        obj.Vertices(i,2) = ymin;
    end
    if abs(obj.Vertices(i,2) - ymax) < 1e-14
        obj.Vertices(i,2) = ymax;
    end
end
ve = [];
for i=1:obj.TotalVertices
    if obj.Vertices(i,1) < xmin || obj.Vertices(i,2) < ymin || obj.Vertices(i,1) > xmax || obj.Vertices(i,2) > ymax
        ve = [ve, i];
    end
end
obj.Vertices(ve,:) = [];
obj.TotalVertices = size(obj.Vertices, 1);
% close file
fclose(fid);
% Allocate Additional Memory Space
n_faces = 0;
face_verts = zeros(0,obj.Dimension);
face_cells = zeros(0,2);
% Generate Cell Vertex Information
for c=1:obj.TotalCells
    c_dg_verts = cell_dg_verts{c};
    dg_verts = vert(c_dg_verts,:);
    obj.CellVerts{c} = zeros(1,length(c_dg_verts));
    % double loop to find cell vertex indices
    for i=1:length(c_dg_verts)
        tv = dg_verts(i,:);
        for j=1:obj.TotalVertices
            if norm(tv - obj.Vertices(j,:)) < 1e-14
                obj.CellVerts{c}(i) = j;
                break
            end
        end
    end
    zind = (obj.CellVerts{c} == 0);
    obj.CellVerts{c}(zind) = [];
    % Determine additional cell information
    obj.CellCenter(c,:) = mean(obj.Vertices(obj.CellVerts{c},:));
    obj.CellVolume(c) = polygonArea(obj.Vertices(obj.CellVerts{c},:));
    sa = 0;
    for i=1:length(obj.CellVerts{c})
        if i==length(obj.CellVerts{c})
            ii = [i,1];
        else
            ii = [i,i+1];
        end
        sa = sa + norm(diff(obj.Vertices(obj.CellVerts{c}(ii),:)));
        obj.CellFaceVerts{c}{i} = obj.CellVerts{c}(ii);
    end
    obj.CellSurfaceArea(c) = sa;
    % Once cell vertices are determine, loop through new vertices and determine
    % faces within the cell. If a new face is discovered, then generate the
    % appropriate information.
    for f = 1:length(obj.CellVerts{c})
        if f==length(obj.CellVerts{c})
            ff = [f,1];
        else
            ff = [f,f+1];
        end
        fverts = obj.CellVerts{c}(ff);
        fverts = sort(fverts);
        ID_v1 = fverts(1);
        ID_v2 = fverts(2);
        ind1 = find( face_verts(:,1) == ID_v1 );
        ind2 = find( face_verts(:,2) == ID_v2 );
        if isempty(ind1) || isempty(ind2)
            n_faces = n_faces + 1;
            face_verts(n_faces,:) = fverts;
            obj.CellFaces{c} = [obj.CellFaces{c},n_faces];
            face_cells(n_faces,:) = [c, 0];
        else
            fbool = true;
            cind = ismember(ind1, ind2);
            for i=1:length(cind)
                if cind(i)
                    fbool = false;
                    fold = ind1(i);
                    break
                end
            end
            if fbool
                n_faces = n_faces + 1;
                face_verts(n_faces,:) = fverts;
                obj.CellFaces{c} = [obj.CellFaces{c},n_faces];
                face_cells(n_faces,:) = [c, 0];
            else
                face_cells(fold,2) = c;
                obj.CellFaces{c} = [obj.CellFaces{c},fold];
            end
        end
    end
end
% Generate Structures
obj.TotalFaces = n_faces;
obj.FaceVerts = cell(obj.TotalFaces,1);
obj.FaceNormal = zeros(obj.TotalFaces,obj.Dimension);
obj.FaceCenter = zeros(obj.TotalFaces,obj.Dimension);
obj.FaceArea = zeros(obj.TotalFaces,1);
obj.FaceCells = uint32(face_cells);
obj.FaceID = zeros(obj.TotalFaces,1,'uint32');
obj.OrthogonalProjection = zeros(obj.TotalFaces, 2);
% Loop through faces and assign face info
for f=1:obj.TotalFaces
    obj.FaceVerts{f} = face_verts(f,:);
    fvv = obj.Vertices(obj.FaceVerts{f},:);
    fc = mean(fvv);
    obj.FaceCenter(f,:) = fc;
    dx = diff(fvv);
    obj.FaceArea(f) = norm(dx);
    obj.FaceNormal(f,:) = [-dx(2),dx(1)];
    obj.FaceNormal(f,:) = obj.FaceNormal(f,:) / norm(obj.FaceNormal(f,:));
    if abs(xmin - fc(1)) < 1e-14 || abs(xmax - fc(1)) < 1e-14 || ...
       abs(ymin - fc(2)) < 1e-14 || abs(ymax - fc(2)) < 1e-14
        obj.FaceID(f) = 1;
    end
    if obj.FaceID(f) == 0
        c1 = face_cells(f,1); c2 = face_cells(f,2);
        obj.OrthogonalProjection(f,1) = get_orthogonal_length(obj.Dimension,obj.CellVolume(c1),length(obj.CellVerts{c1}),obj.FaceArea(f),obj.CellSurfaceArea(c1)); 
        obj.OrthogonalProjection(f,2) = get_orthogonal_length(obj.Dimension,obj.CellVolume(c2),length(obj.CellVerts{c2}),obj.FaceArea(f),obj.CellSurfaceArea(c2));
    else
        c1 = face_cells(f,1);
        obj.OrthogonalProjection(f,1) = get_orthogonal_length(obj.Dimension,obj.CellVolume(c1),length(obj.CellVerts{c1}),obj.FaceArea(f),obj.CellSurfaceArea(c1)); 
    end
    h = obj.OrthogonalProjection(f,1);
    vin = obj.FaceCenter(f,:) + obj.FaceNormal(f,:)*h/1000;
    cv = obj.Vertices(obj.CellVerts{c1},:);
    if inpolygon(vin(1), vin(2), cv(:,1), cv(:,2))
        obj.FaceNormal(f,:) = -1.0*obj.FaceNormal(f,:);
    end
end