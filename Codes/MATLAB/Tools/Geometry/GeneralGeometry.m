%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          General Geometry Mesh Generator
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
classdef GeneralGeometry < handle
    properties (Access = public)
        Dimension
        TotalVertices
        TotalCells
        TotalFaces
        TotalEdges
        TotalInteriorFaces
        TotalBoundaryFaces
        HasPeriodicFaces = false
        IsOrthogonal = false
        IsExtruded = false
    end
    properties (Access = public)
        OriginalMeshType
        MeshType
        Vertices
        
        MatID
        ZoneID
        CellVerts
        CellFaceVerts
        CellCenter
        CellVolume
        CellSurfaceArea
        CellFaces
        CellNeighbors
        CellNeighborFaces
        CellEdges
        
        FaceVerts
        PeriodicFaceVerts
        PeriodicOppositeFaces
        PeriodicFaceCells
        PeriodicBools
        InteriorFaces
        BoundaryFaces
        FaceID
        FaceNormal
        FaceCenter
        FaceArea
        FaceCells
        FaceEdges
        OrthogonalProjection
        
        EdgeVerts
        EdgeCenter
        EdgeLength
        
        VertexCells
        VertexFaces
        CellVertexNumbers
    end
    properties (Access = public) % coord variables
        minX, maxX
        minY, maxY
        minZ, maxZ
        Diameter
    end
    properties (Access = private) % AMR Variables
        OriginalCellCount
        MaxRefinementLevel
        RefinementLevel
        RefinementBool
        RefinementArray
        CellRefinementLevel
        NextCellRefinementLevel
        NumberCellFlags
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Constructor Methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public) % Constructors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GeneralGeometry (varargin)
            n = nargin;
            global glob
            if isempty(glob), glob.print_info = true; end
            if n == 0
                % empty constructor -> do nothing
                return
            elseif n == 1
                if isa(class(varargin{1}), 'CartesianGeometry')
                    obj.transfer_data(varargin{1});
                else
                    error('No clue what you were trying to do here...');
                end
            else
                if glob.print_info, disp('-> Begin General Geometry Construction.'); end
                ttime = tic;
                obj.Dimension = varargin{1};
                % 1D Construction
                if obj.Dimension == 1
                    if ~isa(varargin{2}, 'double')
                        error('1D input needs to consist of a vector of vertex coordinates.')
                    else
                        obj.Constructor_1D(varargin{2});
                        obj.MeshType = 'Quads';
                        obj.OriginalMeshType = 'Quads';
                    end
                    obj.minX = min( obj.Vertices );
                    obj.maxX = max( obj.Vertices );
                    % Calculate Diameter
                    obj.Diameter = obj.maxX - obj.minX;
                % 2D Construction
                elseif obj.Dimension == 2
                    if strcmp(varargin{2}, 'Triangle') || strcmp(varargin{2}, 'triangle')...
                            || strcmp(varargin{2}, 'tri') || strcmp(varargin{2}, 'Tri')
                        [verts, cells, faces] = readTriangleMesh(varargin{3});
                        [tri_out] = process_triangle_file(obj.Dimension,cells,verts,faces);
                        obj.MeshType = 'Triangle';
                        obj.OriginalMeshType = 'Triangle';
                        obj.transfer_data(tri_out);
                    elseif strcmp(varargin{2}, 'Delaunay') || strcmp(varargin{2}, 'delaunay')
                        d_out = process_delaunay_structure_2D(obj.Dimension, varargin{3});
                        obj.MeshType = 'Triangle';
                        obj.OriginalMeshType = 'Triangle';
                        obj.transfer_data(d_out);
                    elseif strcmp(varargin{2}, 'DG') || strcmp(varargin{2}, 'dg')
                        dg_out = process_2D_DG_structure(varargin{3});
                        obj.transfer_data(dg_out);
                        obj.determine_mesh_type()
                        obj.OriginalMeshType = obj.MeshType;
                    end
                    obj.minX = min( obj.Vertices(:,1) );
                    obj.maxX = max( obj.Vertices(:,1) );
                    obj.minY = min( obj.Vertices(:,2) );
                    obj.maxY = max( obj.Vertices(:,2) );
                    % Calculate Diameter
                    dx = obj.maxX - obj.minX;
                    dy = obj.maxY - obj.minY;
                    obj.Diameter = dx; % Default in case of equally-sized domain
                    if dx >= dy
                        obj.Diameter = dx;
                    elseif dy >= dx
                        obj.Diameter = dy;
                    end
                % 3D Construction
                elseif obj.Dimension == 3
                    if strcmp(varargin{2}, 'Tetgen') || strcmp(varargin{2}, 'tetgen')...
                            || strcmp(varargin{2}, 'tet') || strcmp(varargin{2}, 'Tet')
                        [verts, cells, faces, edges] = readTetgenMesh(varargin{3});
                        [tet_out] = process_tetgen_file(obj.Dimension,cells,verts,faces,edges);
                        obj.MeshType = 'Tetrahedron';
                        obj.OriginalMeshType = 'Tetrahedron';
                        obj.transfer_data(tet_out);
                    elseif strcmp(varargin{2}, 'Delaunay') || strcmp(varargin{2}, 'delaunay')
                        d_out = process_delaunay_structure_3D(obj.Dimension, varargin{3});
                        obj.MeshType = 'Tetrahedron';
                        obj.OriginalMeshType = 'Tetrahedron';
                        obj.transfer_data(d_out);
                    end
                    obj.minX = min( obj.Vertices(:,1) );
                    obj.maxX = max( obj.Vertices(:,1) );
                    obj.minY = min( obj.Vertices(:,2) );
                    obj.maxY = max( obj.Vertices(:,2) );
                    obj.minZ = min( obj.Vertices(:,3) );
                    obj.maxZ = max( obj.Vertices(:,3) );
                    % Calculate Diameter
                    dx = obj.maxX - obj.minX;
                    dy = obj.maxY - obj.minY;
                    dz = obj.maxZ - obj.minZ;
                    obj.Diameter = dx; % Default in case of equally-sized domain
                    if dx >= dy && dx >= dz
                        obj.Diameter = dx;
                    elseif dy >= dx && dy >= dz
                        obj.Diameter = dy;
                    elseif dz >= dx && dz >= dy
                        obj.Diameter = dz;
                    end
                end
            end
            
            % Cleanup components
            % ------------------
            obj.determine_faces();
            obj.determine_vertex_components();
            obj.determine_cell_vertex_numbers();
            obj.determine_mesh_type();

            % Populate Initial AMR Variables
            % ------------------------------
            obj.OriginalCellCount = obj.TotalCells;
            obj.NumberCellFlags = 0;
            obj.RefinementBool = 0;
            obj.RefinementLevel = 0;
            obj.MaxRefinementLevel = 0;
            obj.RefinementArray = zeros(obj.TotalCells,2);
            obj.CellRefinementLevel = zeros(obj.TotalCells, 1);

            if glob.print_info
                disp(['-> Total General Generation Time:  ',num2str(toc(ttime))])
                disp(' ')
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Constructor_1D(obj, verts)
            obj.Vertices = verts;
            obj.TotalVertices = length(obj.Vertices);
            obj.TotalCells = obj.TotalVertices - 1;
            obj.TotalFaces = obj.TotalVertices;
            obj.TotalIneriorFaces = obj.TotalVertices - 2;
            obj.TotalBoundaryFaces = 2;
            obj.CellVerts = cell(obj.TotalCells,1);
            obj.CellID = ones(obj.TotalCells,1);
            obj.FaceVerts = cell(obj.TotalCells+1,1);
            for c=1:obj.TotalCells
                obj.CellVerts{c} = [c,c+1];
                obj.CellFaces{c} = [c,c+1];
                obj.FaceVerts{c} = c;
            end
            obj.FaceID = [1,zeros(1,obj.TotalFaces-2),1]';
            obj.FaceVerts{end} = obj.TotalCells + 1;
            obj.FaceArea = ones(length(obj.FaceVerts),1);
            obj.FaceCenter = obj.Vertices;
            obj.FaceNormal = [-1,ones(1,obj.TotalFaces-1)]';
            obj.BoundaryFaces = uint32([1;obj.TotalFaces]);
            obj.InteriorFaces = uint32((2:obj.TotalFaces-1)');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fixup_vertex_orderings(obj)
            if obj.Dimension == 2
                % Loop through cells
                for c=1:obj.TotalCells
                    vind = obj.CellVerts{c};
                    v = obj.Vertices(vind,:);
                    vc = mean(v);
                    [~,ind] = sort(atan2(v(:,2)-vc(2), v(:,1)-vc(1)));
                    nvind = zeros(1,length(vind));
                    for i=1:length(vind)
                        nvind(i) = vind(ind(i));
                    end
                    obj.CellVerts{c} = nvind;
                end
                % Loop through faces
                for f=1:obj.TotalFaces
                    fids = obj.FaceVerts{f};
                    fv = obj.Vertices(fids,:);
                    if obj.FaceID(f) == 0
                        fcells = obj.FaceCells(f,:);
                    else
                        fcells = obj.FaceCells(f,1);
                    end
                    for i=1:length(fcells)
                        vind = obj.CellVerts{fcells(i)};
                        v = obj.Vertices(vind,:);
                        vc = mean(v);
                        [~,ind] = sort(atan2(fv(:,2)-vc(2), fv(:,1)-vc(1)));
                        nvind = zeros(1,length(fids));
                        for j=1:length(fids)
                            nvind(j) = fids(ind(j));
                        end
                        obj.FaceVerts{f} = nvind;
                    end
                end
            elseif obj.Dimension == 3
%                 for f=1:obj.TotalFaces
%                     fids = obj.FaceVerts{f};
%                     fv = obj.Vertices(fids,:);
%                     if obj.FaceID(f) == 0
%                         fcells = obj.FaceCells(f,:);
%                     else
%                         fcells = obj.FaceCells(f,1);
%                     end
%                     
%                 end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_faces(obj)
            IF = 0; BF = 0;
            for f=1:obj.TotalFaces
                if obj.FaceID(f) == 0
                    IF = IF + 1;
                    obj.InteriorFaces(IF) = f;
                else
                    BF = BF + 1;
                    obj.BoundaryFaces(BF) = f;
                end
            end
            obj.InteriorFaces = uint32(obj.InteriorFaces');
            obj.BoundaryFaces = uint32(obj.BoundaryFaces');
            obj.TotalInteriorFaces = length(obj.InteriorFaces);
            obj.TotalBoundaryFaces = length(obj.BoundaryFaces);
            % Determine Face Ordering within Cells
            if obj.Dimension == 2
                for c=1:obj.TotalCells
                    cfaces = obj.CellFaces{c};
                    cverts = obj.CellVerts{c}; nv = length(cverts);
                    tfaces = zeros(1,nv);
                    for i=1:nv
                        if i==nv
                            ii = [i,1];
                        else
                            ii = [i,i+1];
                        end
                        iii = cverts(ii);
                        for ff = 1:nv
                            f = cfaces(ff);
                            fverts = obj.FaceVerts{f};
                            if (iii(1) == fverts(1) && iii(2) == fverts(2)) || (iii(1) == fverts(2) && iii(2) == fverts(1))
                                tfaces(i) = f;
                                break
                            end
                        end
                    end
                    obj.CellFaces{c} = tfaces;
                end
            elseif obj.Dimension == 3
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_vertex_components(obj)
            obj.VertexCells = cell(obj.TotalVertices,1);
            obj.VertexFaces = cell(obj.TotalVertices,1);
            vcbool = ones(obj.TotalVertices,1);
            vfbool = ones(obj.TotalVertices,1);
            vcells = zeros(obj.TotalVertices,2^(obj.Dimension+1));
            vfaces = zeros(obj.TotalVertices,2^(obj.Dimension+1));
            % Loop through cells
            for c=1:obj.TotalCells
                cv = obj.CellVerts{c};
                for i=1:length(cv)
                    ii = cv(i);
                    vcells(ii,vcbool(ii)) = c;
                    vcbool(ii) = vcbool(ii) + 1;
                end
            end
            % Loop through faces
            for f=1:obj.TotalFaces
                fv = obj.FaceVerts{f};
                for i=1:length(fv)
                    ii = fv(i);
                    vfaces(ii,vfbool(ii)) = f;
                    vfbool(ii) = vfbool(ii) + 1;
                end
            end
            % Loop through vertices and assign components
            for v=1:obj.TotalVertices
                obj.VertexCells{v} = vcells(v,1:(vcbool(v)-1));
                obj.VertexFaces{v} = vfaces(v,1:(vfbool(v)-1));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_cell_vertex_numbers(obj)
            obj.CellVertexNumbers = [];
            for c=1:obj.TotalCells
                nv = length(obj.CellVerts{c});
                if length(obj.CellVertexNumbers) < nv
                    obj.CellVertexNumbers(nv) = 1;
                else
                    obj.CellVertexNumbers(nv) = obj.CellVertexNumbers(nv) + 1;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_mesh_type(obj)
            if obj.Dimension == 2
                if length(obj.CellVertexNumbers) == 3
                    obj.MeshType = 'Triangle';
                elseif length(obj.CellVertexNumbers) == 4
                    obj.MeshType = 'Quadrilateral';
                else
                    obj.MeshType = 'Polygon';
                end
            elseif obj.Dimension == 3
                if length(obj.CellVertexNumbers) == 4
                    obj.MeshType = 'Tetrahedron';
                elseif length(obj.CellVertexNumbers) == 8
                    obj.MeshType = 'Hexahedron';
                else
                    obj.MeshType = 'Polyhedron';
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Allocate_Arrays(obj)
            % Cell Arrays
            obj.MatID = ones(obj.TotalCells, 1, 'uint32');
            obj.CellVerts = cell(obj.TotalCells, 1);
            obj.CellFaceVerts = cell(obj.TotalCells, 1);
            obj.CellNeighbors = cell(obj.TotalCells, 1);
            obj.CellNeighborFaces = cell(obj.TotalCells, 1);
            obj.CellCenter = zeros(obj.TotalCells, obj.Dimension);
            obj.CellVolume = zeros(obj.TotalCells, 1);
            obj.CellSurfaceArea = zeros(obj.TotalCells, 1);
            obj.CellFaces = cell(obj.TotalCells, 1);
            % Face Arrays
            obj.FaceVerts = cell(obj.TotalFaces, 1);
            obj.FaceCells = zeros(obj.TotalFaces, 2);
            obj.InteriorFaces = zeros(obj.TotalInteriorFaces,1,'uint32');
            obj.BoundaryFaces = zeros(obj.TotalBoundaryFaces,1,'uint32');
            obj.OrthogonalProjection = zeros(obj.TotalFaces,2);
            obj.FaceID = zeros(obj.TotalFaces, 1, 'uint32');
            obj.FaceNormal = zeros(obj.TotalFaces, obj.Dimension);
            obj.FaceCenter = zeros(obj.TotalFaces, obj.Dimension);
            obj.FaceArea = zeros(obj.TotalFaces, 1);
            % Vertex Arrays
            obj.VertexCells = cell(obj.TotalVertices, 1);
            obj.VertexFaces = cell(obj.TotalVertices, 1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function allocate_more_memory(obj, nverts, ncells, nfaces)
            % Cell Arrays
            obj.MatID = [obj.MatID;zeros(ncells,1,'uint32')];
            obj.CellVerts = [obj.CellVerts;cell(ncells,1)];
            obj.CellNeighbors = [obj.CellNeighbors;cell(ncells,1)];
            obj.CellNeighborFaces = [obj.CellNeighborFaces;cell(ncells,1)];
            obj.CellCenter = [obj.CellCenter;zeros(ncells,obj.Dimension)];
            obj.CellVolume = [obj.CellVolume;zeros(ncells,1)];
            obj.CellSurfaceArea = [obj.CellSurfaceArea;zeros(ncells,1)];
            obj.CellFaces = [obj.CellFaces;cell(ncells,1)];
            % Face Arrays
            obj.FaceVerts = [obj.FaceVerts;cell(nfaces,1)];
            obj.FaceCells = [obj.FaceCells;zeros(nfaces,2)];
            obj.OrthogonalProjection = [obj.OrthogonalProjection;zeros(nfaces,2)];
            obj.FaceID = [obj.FaceID;zeros(nfaces,1)];
            obj.FaceNormal = [obj.FaceNormal;zeros(nfaces,obj.Dimension)];
            obj.FaceCenter = [obj.FaceCenter;zeros(nfaces,obj.Dimension)];
            obj.FaceArea = [obj.FaceArea;zeros(nfaces,1)];
            % Vertex Arrays
            obj.VertexCells = [obj.VertexCells;cell(nverts,1)];
            obj.VertexFaces = [obj.VertexFaces;cell(nverts,1)];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calculate_orthogonal_projections(obj)
            for c=1:obj.TotalCells
                cf = obj.CellFaces{c}; ncf = length(cf);
                cv = obj.CellVerts{c}; ncv = length(cv);
                fcnodes = cell(ncf, 1);
                for f=1:ncf
                    ff = cf(f);
                    fverts = obj.FaceVerts{ff}; nfverts = length(fverts);
                    tfn = zeros(1, nfverts);
                    for i=1:nfverts
                        for j=1:ncv
                            if fverts(i) == cv(j);
                                tfn(i) = j;
                                break
                            end
                        end
                    end
                    fcnodes{f} = tfn;
                end
                v = obj.Vertices(cv,:);
                h = get_orthogonal_projection(v,fcnodes);
                % Place projection into array
                for f=1:ncf
                    ff = cf(f);
                    fcells = obj.FaceCells(ff,:);
                    if fcells(1) == c
                        obj.OrthogonalProjection(ff,1) = h(f);
                    else
                        obj.OrthogonalProjection(ff,2) = h(f);
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Accessor Routines - analogous to 'get' commands in C++
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_vertices(obj)
            out = obj.Vertices;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_faces(obj)
            out = obj.FaceVerts;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_cell_verts(obj)
            out = obj.CellVerts;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_all_cell_faces(obj)
            out = obj.CellFaces;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_material_id(obj, cellID)
            out = obj.MatID(cellID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_verts(obj, cellID)
            out = obj.Vertices(obj.CellVerts{cellID},:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_vert_indices(obj, cellID)
            out = obj.CellVerts{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_faces(obj, cellID)
            out = obj.CellFaces{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_centers(obj, cellID)
            out = obj.CellCenter(cellID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_volumes(obj, cellID)
            out = obj.CellVolumes(cellID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_cell_surface_areas(obj, cellID)
            out = obj.CellSurfaceArea(cellID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_interior_face(obj, x)
            out = obj.InteriorFaces(x);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_boundary_face(obj, x)
            out = obj.BoundaryFaces(x);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_verts(obj, faceID)
            out = obj.Vertices(obj.FaceVerts{faceID},:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_vert_indices(obj, faceID)
            out = obj.FaceVerts{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_cells(obj, faceID)
            out = obj.FaceCells(faceID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_orth_proj(obj, faceID, cellID)
            out = obj.OrthogonalProjection(faceID,cellID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_normals(obj, faceID)
            out = obj.FaceNormal(faceID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_areas(obj, faceID)
            out = obj.FaceArea(faceID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_flags(obj, faceID)
            out = obj.FaceID(faceID);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_face_centers(obj, faceID)
            out = obj.FaceCenter(faceID,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_mesh_type(obj)
            out = obj.MeshType;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_refinement_level(obj)
            out = obj.RefinementLevel;
        end
    end
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                   List of Extended Accessor Routines
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_triangular_connectivity_array(obj)
            if obj.Dimension == 2
                if ~strcmp(obj.MeshType, 'Triangle')
                    error('Can only get 2D connectivity array for triangles.')
                end
            elseif obj.Dimension == 3
                if ~strcmp(obj.MeshType, 'Tetrahedron')
                    error('Can only get 3D connectivity array for tetrahedra.')
                end
            end
            out = zeros(obj.TotalCells, obj.Dimension+1);
            for c=1:obj.TotalCells
                out(c,:) = obj.CellVerts{c};
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Generation Routines (i.e. extrude, set flags, randomize)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_cell_matIDs_inside_domain(obj, val, verts, faces)
            % Function to loop through cells and determine if cell
            % center lies within shape defined by 'verts' and 'faces'.
            % 'faces' is only needed in 3 dimensions.
            % -------------------------------------------------------------
            if obj.Dimension == 1
                for c=1:obj.TotalCells
                    if obj.CellCenter(c) > verts(1) && obj.CellCenter(c) < verts(2)
                        obj.MatID(c) = val;
                    end
                end
            elseif obj.Dimension == 2
                [in] = inpoly(obj.CellCenter,verts);
                for c=1:obj.TotalCells
                    if logical(in(c))
                        obj.MatID(c) = val;
                    end
                end
            elseif obj.Dimension == 3
                FV.verts = verts; FV.faces = faces;
                bools = inpolyhedron(FV,obj.CellCenter);
                for c=1:obj.TotalCells
                    if logical(bools(c))
                        obj.MatID(c) = val;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_face_flag_on_surface(obj, val, verts)
            if obj.Dimension == 1
                for f=1:obj.TotalFaces
                    cf = obj.FaceCenter(f);
                    if abs(cf - verts) < 1e-14
                        obj.FaceID(f) = val;
                    end
                end
            elseif obj.Dimension == 2
                for f=1:obj.TotalFaces
                    cf = obj.FaceCenter(f,:);
                    v1 = verts(1,:)-cf;
                    v2 = verts(2,:)-cf;
                    vv = verts(2,:) - verts(1,:);
                    mat = [v1;v2];
                    if abs(det(mat)) < 1e-14
                        Kp = (vv)*(-v1)';
                        K12 = vv*vv';
                        if Kp > 0 && Kp < K12
                            obj.FaceID(f) = val;
                        end
                    end
                end
            elseif obj.Dimension == 3
                nfv = size(verts,1);
                for f=1:obj.TotalFaces
                    fc = obj.FaceCenter(f,:);
                    phi_sum = 0;
                    for i=1:nfv
                        ii = [i,mod(i,nfv)+1];
                        vt = [verts(ii,:);fc];
                        p1 = vt(1,:) - fc;
                        p2 = vt(2,:) - fc;
                        n1 = norm(p1);
                        n2 = norm(p2);
                        cphi = (p1*p2') / (n1*n2);
                        phi_sum = phi_sum + acos(cphi);
                    end
                    if abs(2*pi - phi_sum) < 1e-14
                        obj.FaceID(f) = val;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_all_boundary_flags(obj, val)
            obj.FaceID(obj.BoundaryFaces) = val;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_periodic_flag(obj, val, dim)
            global glob
            obj.HasPeriodicFaces = true;
            if isempty(obj.PeriodicFaceVerts) || length(obj.PeriodicFaceVerts) ~= obj.TotalFaces
                obj.PeriodicFaceVerts = cell(obj.TotalFaces, 1); 
            end
            if isempty(obj.PeriodicBools) || length(obj.PeriodicFaceVerts) ~= obj.TotalFaces
                obj.PeriodicBools = logical(zeros(obj.TotalFaces, 1)); 
            end
            if isempty(obj.PeriodicFaceCells) || length(obj.PeriodicFaceCells) ~= obj.TotalFaces
                obj.PeriodicFaceCells = zeros(obj.TotalFaces, 1); 
            end
            if isempty(obj.PeriodicOppositeFaces) || length(obj.PeriodicOppositeFaces) ~= obj.TotalFaces
                obj.PeriodicOppositeFaces = zeros(obj.TotalFaces, 1); 
            end
            if obj.Dimension == 1
                obj.FaceID(1) = val;
                obj.FaceID(end) = val;
                obj.PeriodicFaceVerts{1} = [1,obj.TotalVertices];
                obj.PeriodicFaceVerts{obj.TotalVertices} = [1,1];
                obj.PeriodicOppositeFaces(1) = obj.TotalVertices;
                obj.PeriodicOppositeFaces(end) = 1;
                obj.PeriodicFaceCells(1) = obj.TotalCells;
                obj.PeriodicFaceCells(end) = 1;
                obj.PeriodicBools([1, end]) = true;
                return
            else
                if ~isa(dim, 'char'), error('Periodic Conditions requires either x, y, or z.'); end
                dim = lower(dim);
                if strcmp( dim, 'x' )
                    ind = 1; dbnds = [min(obj.Vertices(:,1)), max(obj.Vertices(:,1))];
                elseif strcmp( dim, 'y' )
                    if obj.Dimension < 1, error('No y direction in 1D.'); end
                    ind = 2; dbnds = [min(obj.Vertices(:,2)), max(obj.Vertices(:,2))];
                elseif strcmp( dim, 'z' )
                    if obj.Dimension ~= 3, error('No z direction in 1D/2D.'); end
                    ind = 3; dbnds = [min(obj.Vertices(:,3)), max(obj.Vertices(:,3))];
                end
                oind = 1:obj.Dimension; oind(ind) = [];
            end
            for f=1:obj.TotalBoundaryFaces
                bf = obj.BoundaryFaces(f);
                if obj.PeriodicBools(bf), continue; end
                % Determine if face is on minimum dimension
                if abs( obj.FaceCenter(bf,ind) - dbnds(1) ) < glob.small
                    fb_cent = obj.FaceCenter(bf,oind);
                    % Find matching face on maximum dimension
                    tbool = false;
                    for ff=1:obj.TotalBoundaryFaces
                        bff = obj.BoundaryFaces(ff);
                        if abs( obj.FaceCenter(bff,ind) - dbnds(2) ) < glob.small && ...
                           norm( obj.FaceCenter(bff,oind) - fb_cent) < glob.small
                            tbool = true;
                            break;
                        end
                    end
                    if ~tbool, error('Cell faces do not align for periodic boundary conditions.'); end
                    obj.FaceID([bf, bff]) = val;
                    obj.PeriodicOppositeFaces(bf) = bff;
                    obj.PeriodicOppositeFaces(bff) = bf;
                    obj.PeriodicFaceCells(bf) = obj.FaceCells(bff,1);
                    obj.PeriodicFaceCells(bff) = obj.FaceCells(bf,1);
                    obj.PeriodicBools([bf, bff]) = true;
                    vfind = obj.FaceVerts{bf}; nv = length(vfind);
                    vffind = obj.FaceVerts{bff};
                    obj.PeriodicFaceVerts{bf} = zeros(nv,2);
                    obj.PeriodicFaceVerts{bff} = zeros(nv,2);
                    fverts = obj.Vertices(vfind,:);
                    ffverts = obj.Vertices(vffind,:);
                    for i=1:nv
                        for j=1:nv
                            if norm(fverts(i,oind) - ffverts(j,oind)) < glob.small
                                obj.PeriodicFaceVerts{bf}(i,:) = [j, vffind(j)];
                                obj.PeriodicFaceVerts{bff}(j,:) = [i, vfind(i)];
                                break
                            end
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function remove_random_verices(obj, N)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function extrude_mesh_2D_to_3D(obj, levels)
            ttime = tic;
            global glob
            if glob.print_info, disp(['-> Begin Geometry Extrusion from ',num2str(obj.Dimension),'D to ',num2str(obj.Dimension+1),'D']); end
            extrude_mesh(obj,levels);
            % Cleanup Mesh Information
            obj.determine_faces();
            obj.determine_vertex_components();
            obj.determine_cell_vertex_numbers();
            obj.determine_mesh_type();
            if glob.print_info, disp(['-> Total Geometry Extrusion Time:  ',num2str(toc(ttime))]); end
            if glob.print_info, disp(' '); end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Unstructured Routines - this is not tested for AMR...
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function clear_refinement_array(obj)
            obj.RefinementArray = zeros(obj.TotalCells,2,'uint32');
            obj.RefinementBool = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_refinement_flag(obj, cellID, flag, val)
            obj.RefinementBool = 1;
            if nargin < 3
                obj.RefinementArray(cellID,:) = [flag, 1];
            else
                obj.RefinementArray(cellID,:) = [flag, val];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function refine_mesh(obj)
            for c=1:obj.TotalCells
                if obj.RefinementArray(c,1) ~= 0
                    refine_cell(c);
                end
            end
            obj.RefinementLevel = obj.RefinementLevel + 1;
            obj.clear_refinement_array();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Private Functions for Internal Calculations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function transfer_data(obj, obj_in)
            obj.TotalCells = obj_in.TotalCells;
            obj.TotalFaces = obj_in.TotalFaces;
            obj.Vertices = obj_in.Vertices;
            obj.TotalVertices = size(obj.Vertices,1);
            obj.CellVerts = obj_in.CellVerts;
            obj.CellFaces = obj_in.CellFaces;
            obj.CellFaceVerts = obj_in.CellFaceVerts;
            obj.CellCenter = obj_in.CellCenter;
            obj.CellVolume = obj_in.CellVolume;
            obj.CellSurfaceArea = obj_in.CellSurfaceArea;
            obj.FaceVerts = obj_in.FaceVerts;
            obj.MatID = obj_in.MatID;
            obj.FaceID = obj_in.FaceID;
            obj.OrthogonalProjection = obj_in.OrthogonalProjection;
            obj.FaceCenter = obj_in.FaceCenter;
            obj.FaceNormal = obj_in.FaceNormal;
            obj.FaceArea = obj_in.FaceArea;
            obj.FaceCells = obj_in.FaceCells;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           Accessory Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_cell(obj, cellID)
if obj.Dimension == 1
    refine_cell_1D(obj, cellID);
elseif obj.Dimension == 2
    refine_cell_2D(obj, cellID);
elseif obj.Dimension == 3
    refine_cell_3D(obj, cellID);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_cell_1D(obj, cellID)
val = obj.RefinementArray(cellID,2);
cverts = obj.CellVerts{cellID};
verts = obj.Vertices(cverts,:);
dx = abs(verts(1) - verts(2));
nc = 2^val;
ddx = dx/nc;
newc = nc-1;
for r=1:nc
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_cell_2D(obj, cellID)
flag = obj.RefinementArray(cellID,1);
val = obj.RefinementArray(cellID,2);
switch(flag)
    case(1) % Add new point in cell center
        verts = obj.Vertices(obj.CellVerts{cellID},:);
        rcenter = mean(verts);
        faces = obj.CellFaces{cellID};
        newcells = length(faces)-1;
        newfaces = newcells;
        % Loop through nverts-1
        for i=1:newcells
            nverts = [verts([i,i+1],:);rcenter];
            
        end
        % Last cell addition
        
        obj.TotalCells = obj.TotalCells + newcells;
    case(2)
        verts = obj.Vertices(obj.CellVerts{cellID},:);
        if size(verts,1) == 3
            ncells = 4;
        elseif size(verts,1) == 4
            ncells = 4;
        else
            error('This refinement option only works for triangles and quads.')
        end
    otherwise
        
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_cell_3D(obj, cellID)
flag = obj.RefinementArray(cellID,1);
val = obj.RefinementArray(cellID,2);
switch(flag)
    case(1)
        
    case(2)
        
    otherwise
        
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function extrude_mesh(obj, levels)
global glob
% Get new Geometry Information
% ----------------------------
new_dim = obj.Dimension + 1;
nlev = length(levels);
num_new_verts = obj.TotalVertices * nlev;
num_new_cells = obj.TotalCells * (nlev - 1);
num_new_faces = obj.TotalCells * nlev + obj.TotalFaces * (nlev - 1);
v_stride  = obj.TotalVertices;
c_stride  = obj.TotalCells;
fe_stride = obj.TotalFaces;
fb_stride = obj.TotalCells;
ef_stride = obj.TotalFaces;
ev_stride = obj.TotalVertices;
ef_tot = ef_stride*nlev;
ev_tot = ev_stride*(nlev-1);
num_ef = obj.TotalFaces*nlev;
num_ev = obj.TotalVertices*(nlev-1);
num_fe = obj.TotalFaces * (nlev - 1);
num_fb = obj.TotalCells * nlev;
% Allocate New Memory Information
% -------------------------------
new_verts = zeros(num_new_verts, new_dim);
new_cell_verts = cell(num_new_cells, 1);
new_face_verts = cell(num_new_faces, 1);
new_cell_faces = cell(num_new_cells, 1);
new_face_cells = zeros(num_new_faces, 2);
new_cell_matids = zeros(num_new_cells, 1);
new_cell_vols = zeros(num_new_cells, 1);
new_cell_sa = zeros(num_new_cells, 1);
new_cell_centers = zeros(num_new_cells, new_dim);
new_face_ids = zeros(num_new_faces, 1);
new_face_areas = zeros(num_new_faces, 1);
new_face_centers = zeros(num_new_faces, new_dim);
new_face_norms = zeros(num_new_faces, new_dim);
new_orth_len = zeros(num_new_faces, 2);
if new_dim == 3
    num_new_edges = obj.TotalVertices*(nlev-1) + obj.TotalFaces*nlev;
    new_edge_verts = zeros(num_new_edges,2);
    new_edge_centers = zeros(num_new_edges,new_dim);
    new_edge_lengths = zeros(num_new_edges,1);
    new_cell_edges = cell(num_new_cells, 1);
    new_face_edges = cell(num_new_faces, 1);
end
% Loop through levels and get new vertex info
for i=1:nlev
    vv = (i-1)*v_stride+1:i*v_stride;
    new_verts(vv,:) = [obj.Vertices,levels(i)*ones(obj.TotalVertices, 1)];
    if new_dim == 3
        ef = (i-1)*ef_stride+1:i*ef_stride;
        for f=1:obj.TotalFaces
            new_edge_verts((i-1)*obj.TotalFaces+f,:) = vv(obj.FaceVerts{f});
        end
        new_edge_centers(ef,:) = [obj.FaceCenter,levels(i)*ones(obj.TotalFaces,1)];
        new_edge_lengths(ef) = obj.FaceArea;
        if i ~=nlev
            ev = num_ef+(i-1)*ev_stride+1:num_ef+i*ev_stride;
            vv2 = i*v_stride+1:(i+1)*v_stride;
            new_edge_verts(ev,:) = [vv',vv2'];
            new_edge_centers(ev,:) = [obj.Vertices,(levels(i+1) + levels(i))/2*ones(obj.TotalVertices, 1)];
            new_edge_lengths(ev) = levels(i+1) - levels(i);
        end
    end
end
% Loop through levels and get new cell/face info
if glob.print_info, disp('   -> Begin New Cell and Face Generation.'); end
cc = 0; ff = 0;
for i=1:nlev-1
    vv1 = (i-1)*v_stride+1:i*v_stride;
    vv2 = i*v_stride+1:(i+1)*v_stride;
    ffe = (i-1)*fe_stride+1:i*fe_stride;
    cclv = (i-1)*c_stride+1:i*c_stride;
    fb1 = (i-1)*fb_stride + num_fe + 1:i*fb_stride + num_fe;
    fb2 = i*fb_stride + num_fe + 1:(i+1)*fb_stride + num_fe;
    n_center = mean(levels([i,i+1]));
    if new_dim == 3
        ef1 = (i-1)*ef_stride+1:i*ef_stride;
        ef2 = i*ef_stride+1:(i+1)*ef_stride;
        ev = num_ef+(i-1)*ev_stride+1:num_ef+i*ev_stride;
    end
    % Loop through faces
    new_face_ids(ffe) = obj.FaceID;
    new_face_centers(ffe,:) = [obj.FaceCenter, ones(fe_stride,1)*n_center];
    new_face_norms(ffe,:) = [obj.FaceNormal, zeros(fe_stride,1)];
    for f=1:obj.TotalFaces
        fverts = obj.FaceVerts{f};
        cverts = obj.CellVerts{obj.FaceCells(f,1)};
        for j=1:length(cverts)
            if fverts(1) == cverts(j)
                if j==length(cverts)
                    if fverts(2) == cverts(1)
                        sbool = true;
                    else
                        sbool = false;
                    end
                else
                    if fverts(2) == cverts(j+1);
                        sbool = true;
                    else
                        sbool = false;
                    end
                end
                break
            end
        end
        ff = ffe(f);
        if sbool
            new_face_verts{ff} = [fliplr(vv1(fverts)), (vv2(fverts))];
        else
            new_face_verts{ff} = [vv1(fverts), fliplr(vv2(fverts))];
        end
        if new_dim == 3
            new_face_areas(ff) = polygonArea3d(new_verts(new_face_verts{ff},:));
        elseif new_dim == 2
            new_face_areas(ff) = norm(diff(new_verts(new_face_verts{ff},:)));
        end
        if new_face_ids(ff) == 0 % interior
            new_face_cells(ff,:) = cclv(obj.FaceCells(f,:));
        else % boundary
            new_face_cells(ff,1) = cclv(obj.FaceCells(f,1));
        end
        % Assign 3D edges to face
        if new_dim == 3
            new_face_edges{ff} = [ef1(f),ev(fverts(1)),ef2(f),ev(fverts(2))];
        end
    end
    % Loop through cells
    for c=1:obj.TotalCells
        cc = cc + 1;
        ccvs = obj.CellVerts{c};
        cfaces = obj.CellFaces{c};
        new_cell_verts{cc} = [vv1(ccvs),vv2(ccvs)];
        new_cell_matids(cc) = obj.MatID(c);
        new_cell_centers(cc,:) = mean(new_verts(new_cell_verts{cc}, :));
        new_cell_faces{cc} = [ffe(cfaces), fb1(c), fb2(c)];
        new_cell_sa(cc) = sum(new_face_areas(ffe(cfaces))) + 2*obj.CellVolume(c);
        new_cell_vols(cc) = obj.CellVolume(c) * (levels(i+1) - levels(i));
        % Assign 3D edges to cell
        if new_dim == 3
            new_cell_edges{cc} = [ef1(cfaces), ef2(cfaces), ev(ccvs)];
        end
    end
end
% Loop through all levels
if glob.print_info, disp('   -> Begin Cell and Face Cleanup Operations.'); end
for i=1:nlev
    fb = num_fe+(i-1)*fb_stride+1:num_fe+i*fb_stride;
    vv = (i-1)*v_stride+1:i*v_stride;
    cclv1 = (i-2)*c_stride+1:(i-1)*c_stride;
    cclv2 = (i-1)*c_stride+1:i*c_stride;
    new_face_centers(fb,:) = [obj.CellCenter, ones(length(fb),1)*levels(i)];
    new_face_areas(fb) = obj.CellVolume;
    % Assign 3D edges to faces
    if new_dim == 3
        ef = (i-1)*ef_stride+1:i*ef_stride;
        for c=1:obj.TotalCells
            new_face_edges{fb(c)} = ef(obj.CellFaces{c});
        end
    end
    for c=1:obj.TotalCells
        ccvs = obj.CellVerts{c};
        if i==1
            new_face_verts{fb(c)} = vv(ccvs);
        else
            new_face_verts{fb(c)} = fliplr(vv(ccvs));
        end
    end
    if i==1
        new_face_norms(fb,:) = [zeros(length(fb),obj.Dimension), -1*ones(length(fb),1)];
        new_face_ids(fb) = 1;
    elseif i==nlev
        new_face_norms(fb,:) = [zeros(length(fb),obj.Dimension), ones(length(fb),1)];
        new_face_ids(fb) = 1;
    else
        new_face_norms(fb,:) = [zeros(length(fb),obj.Dimension), ones(length(fb),1)];
        new_face_ids(fb) = 0;
    end
    for c=1:obj.TotalCells
        ffb = fb(c);
        if i==1
            new_face_cells(ffb,1) = cclv2(c);
        elseif i==nlev
            new_face_cells(ffb,1) = cclv1(c);
        else
            new_face_cells(ffb,:) = [cclv1(c), cclv2(c)];
        end
    end
end
% Loop through faces for orthogonal projections
if glob.print_info, disp('   -> Begin Face Orthogonal Projection Calculations.'); end
for f=1:num_new_faces
    fcells = new_face_cells(f,:);
    if fcells(2) == 0, fcells(2) = []; end
    for i=1:length(fcells)
        nv = length(new_cell_verts{fcells(i)});
        new_orth_len(f,i) = get_orthogonal_length(new_dim, new_cell_vols(fcells(i)), nv, new_face_areas(f), new_cell_sa(fcells(i)));
    end
end
% Set New Information
obj.IsExtruded = true;
obj.Dimension = new_dim;
obj.TotalVertices = num_new_verts;
obj.TotalCells = num_new_cells;
obj.TotalFaces = num_new_faces;
obj.Vertices = new_verts;
obj.CellVerts = new_cell_verts;
obj.CellFaces = new_cell_faces;
obj.MatID = new_cell_matids;
obj.CellVolume = new_cell_vols;
obj.CellCenter = new_cell_centers;
obj.CellSurfaceArea = new_cell_sa;
obj.OrthogonalProjection = new_orth_len;
obj.FaceVerts = new_face_verts;
obj.FaceCells = new_face_cells;
obj.FaceID = new_face_ids;
obj.FaceArea = new_face_areas;
obj.FaceNormal = new_face_norms;
obj.FaceCenter = new_face_centers;
if new_dim == 3
    obj.EdgeVerts = new_edge_verts;
    obj.EdgeCenter = new_edge_centers;
    obj.EdgeLength = new_edge_lengths;
    obj.CellEdges = new_cell_edges;
    obj.FaceEdges = new_face_edges;
end
end

