%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          AMR Geometry Mesh Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB class to generate all data structures necessary
%                   to fully describe a geometric domain to be used for
%                   finite element (FEM) calculations.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef AMRGeometry < handle
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
        CellNeighbors
        CellNeighborFaces
        CellCenter
        CellVolume
        CellSurfaceArea
        CellFaces
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
    properties (Access = public) % AMR variables
        MaxIrregularity = inf
        OriginalCellCount
        OriginalFaceCount
        OriginalVertexCount
        NewCellsPerRefinement
        OriginalFacesPerCell
        MeshRefinementLevel
        PreviousCell
        CellRefinedLastCycle
        CellRefinementFlag
        CellRefinementLevel
        CellRefinementTree
        CellRefinementTop
        CellRefinementTreeHierarchy
        CellRefinementTreeHierarchyLevel
        AdjTreeCells
        AdjTreeCellFaces
        AdjTreeCellFaceNumbering
        AdjTreeFaceCells
        NextCellRefinementLevel
        RefCellCornerVerts
        RefCellFaceVerts
        RefCellMidFaceVerts
        RefCellFaces
        RefCellFaceCells
        RefCellFaceNumbering
        RefCellHigherLvls
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Constructor Methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = AMRGeometry (varargin)
            n = nargin;
            global glob
            if isempty(glob), glob.print_info = true; end
            if n == 0
                % empty constructor -> do nothing
                return
            else
                
                % Quick error checking
                if n > 1, error('Input is simply a prior Geometry object.'); end
                if      ~isa(varargin{1},'GeneralGeometry') && ...
                        ~isa(varargin{1},'CartesianGeometry')    
                    error('Input requires prior Geometry object.')
                end
                if      strcmp('Polygon',varargin{1}.OriginalMeshType) || ...
                        strcmp('Polyhedron',varargin{1}.OriginalMeshType)
                    error('Currently cannot perform AMR on polytope meshes.')
                end
                % Prepare Input Space
                mesh = varargin{1}; clear varargin;
                obj.Dimension               = mesh.Dimension;
                obj.TotalVertices           = mesh.TotalVertices;
                obj.TotalCells              = mesh.TotalCells;
                obj.TotalFaces              = mesh.TotalFaces;
                obj.TotalEdges              = mesh.TotalEdges;
                obj.TotalInteriorFaces      = mesh.TotalInteriorFaces;
                obj.TotalBoundaryFaces      = mesh.TotalBoundaryFaces;
                obj.HasPeriodicFaces        = mesh.HasPeriodicFaces;
                obj.IsOrthogonal            = mesh.IsOrthogonal;
                obj.IsExtruded              = mesh.IsExtruded;
                obj.OriginalMeshType        = mesh.OriginalMeshType;
                obj.MeshType                = mesh.MeshType;
                obj.Vertices                = mesh.Vertices;
                obj.MatID                   = mesh.MatID;
                obj.ZoneID                  = mesh.ZoneID;
                obj.CellVerts               = mesh.CellVerts;
                obj.CellFaceVerts           = mesh.CellFaceVerts;
                obj.CellNeighbors           = mesh.CellNeighbors;
                obj.CellNeighborFaces       = mesh.CellNeighborFaces;
                obj.CellCenter              = mesh.CellCenter;
                obj.CellVolume              = mesh.CellVolume;
                obj.CellSurfaceArea         = mesh.CellSurfaceArea;
                obj.CellFaces               = mesh.CellFaces;
                obj.CellEdges               = mesh.CellEdges;
                obj.FaceVerts               = mesh.FaceVerts;
                obj.PeriodicFaceVerts       = mesh.PeriodicFaceVerts;
                obj.PeriodicOppositeFaces   = mesh.PeriodicOppositeFaces;
                obj.PeriodicFaceCells       = mesh.PeriodicFaceCells;
                obj.PeriodicBools           = mesh.PeriodicBools;
                obj.InteriorFaces           = mesh.InteriorFaces;
                obj.BoundaryFaces           = mesh.BoundaryFaces;
                obj.FaceID                  = mesh.FaceID;
                obj.FaceNormal              = mesh.FaceNormal;
                obj.FaceCenter              = mesh.FaceCenter;
                obj.FaceArea                = mesh.FaceArea;
                obj.FaceCells               = mesh.FaceCells;
                obj.FaceEdges               = mesh.FaceEdges;
                obj.OrthogonalProjection    = mesh.OrthogonalProjection;
                obj.EdgeVerts               = mesh.EdgeVerts;
                obj.EdgeCenter              = mesh.EdgeCenter;
                obj.EdgeLength              = mesh.EdgeLength;
                obj.VertexCells             = mesh.VertexCells;
                obj.VertexFaces             = mesh.VertexFaces;
                obj.CellVertexNumbers       = mesh.CellVertexNumbers;
                obj.minX                    = mesh.minX;
                obj.maxX                    = mesh.maxX;
                obj.minY                    = mesh.minY;
                obj.maxY                    = mesh.maxY;
                obj.minZ                    = mesh.minZ;
                obj.maxZ                    = mesh.maxZ;
                obj.Diameter                = mesh.Diameter;
                % Other Preliminary Information
                if isa(mesh, 'CartesianGeometry')
                    obj.NewCellsPerRefinement = 2^obj.Dimension-1;
                    obj.OriginalFacesPerCell = 2^obj.Dimension;
                elseif isa(mesh, 'GeneralGeometry')
                    obj.NewCellsPerRefinement = obj.Dimension + 1;
                    obj.OriginalFacesPerCell = obj.Dimension + 1;
                end
                obj.OriginalCellCount = obj.TotalCells;
                obj.OriginalFaceCount = obj.TotalFaces;
                obj.OriginalVertexCount = obj.TotalVertices;
                obj.MeshRefinementLevel = 0;
                obj.PreviousCell = (1:obj.TotalCells)';
                obj.CellRefinementFlag = false(obj.TotalCells, 1);
                obj.CellRefinementLevel = zeros(obj.TotalCells, 1);
                obj.CellRefinementTree = cell(obj.TotalCells, 1);
                obj.CellRefinementTop = (1:obj.TotalCells)';
                obj.CellRefinementTreeHierarchy = cell(obj.TotalCells, 1);
                obj.CellRefinementTreeHierarchyLevel = zeros(obj.TotalCells,1);
                obj.AdjTreeCells = cell(obj.TotalCells, 1);
                obj.AdjTreeCellFaceNumbering = cell(obj.TotalCells, 1);
                obj.AdjTreeCellFaces = obj.CellFaces;
                obj.AdjTreeFaceCells = obj.FaceCells;
                obj.RefCellCornerVerts = obj.CellVerts;
                obj.RefCellFaceVerts = obj.CellFaceVerts;
                obj.RefCellMidFaceVerts = cell(obj.TotalCells, 1);
                obj.RefCellFaces = cell(obj.TotalCells, 1);
                obj.RefCellFaceCells = cell(obj.TotalCells, 1);
                obj.RefCellFaceNumbering = cell(obj.TotalCells, 1);
                obj.RefCellHigherLvls = cell(obj.TotalCells, 1);
                for c=1:obj.TotalCells
                    obj.CellRefinementTree{c} = c;
                    obj.CellRefinementTreeHierarchy{c} = c;
                    obj.AdjTreeCells{c} = zeros(1,obj.OriginalFacesPerCell);
                    obj.AdjTreeCellFaceNumbering{c} = zeros(1,obj.OriginalFacesPerCell);
                    cfaces = obj.CellFaces{c};
                    obj.RefCellFaces{c} = cell(length(cfaces),1);
                    obj.RefCellFaceCells{c} = cell(length(cfaces),1);
                    obj.RefCellMidFaceVerts{c} = zeros(1,length(cfaces));
                    obj.RefCellFaceNumbering{c} = cell(1,length(cfaces));
                    obj.RefCellHigherLvls{c} = false(1,length(cfaces));
                    for ff=1:length(cfaces)
                        f = cfaces(ff);
                        obj.RefCellFaces{c}{ff} = f;
                        fid = obj.FaceID(f);
                        if fid == 0
                            fcells = obj.FaceCells(f,:);
                            if c==fcells(1)
                                obj.AdjTreeCells{c}(ff) = fcells(2);
                                obj.RefCellFaceCells{c}{ff} = fcells(2);
                            elseif c==fcells(2)
                                obj.AdjTreeCells{c}(ff) = fcells(1);
                                obj.RefCellFaceCells{c}{ff} = fcells(1);
                            end
                            cfaces2 = obj.CellFaces{obj.AdjTreeCells{c}(ff)};
                            for ff2=1:length(cfaces2)
                                if f==cfaces2(ff2)
                                    obj.AdjTreeCellFaceNumbering{c}(ff) = ff2;
                                    obj.RefCellFaceNumbering{c}{ff} = ff2;
                                end
                            end
                        end
                    end
                end
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
        function allocate_more_memory(obj, nverts, ncells, nfaces)
            % Cell Arrays
            if ncells > 0
                new_cells = cell(ncells,1);
                obj.MatID = [obj.MatID;zeros(ncells,1,'uint32')];
                obj.CellVerts = [obj.CellVerts;cell(ncells,1)];
                obj.CellFaceVerts = [obj.CellFaceVerts;new_cells];
                obj.CellNeighbors = [obj.CellNeighbors;new_cells];
                obj.CellNeighborFaces = [obj.CellNeighborFaces;new_cells];
                obj.CellCenter = [obj.CellCenter;zeros(ncells,obj.Dimension)];
                obj.CellVolume = [obj.CellVolume;zeros(ncells,1)];
                obj.CellSurfaceArea = [obj.CellSurfaceArea;zeros(ncells,1)];
                obj.CellFaces = [obj.CellFaces;new_cells];
                obj.PreviousCell = [obj.PreviousCell;zeros(ncells,1)];
                obj.CellRefinementFlag = [obj.CellRefinementFlag;false(ncells,1)];
                obj.CellRefinementLevel = [obj.CellRefinementLevel;zeros(ncells,1)];
                obj.CellRefinementTop = [obj.CellRefinementTop;zeros(ncells,1)];
                obj.CellRefinementTreeHierarchy = [obj.CellRefinementTreeHierarchy;new_cells];
                obj.RefCellFaces = [obj.RefCellFaces;new_cells];
                obj.RefCellFaceCells = [obj.RefCellFaceCells;new_cells];
                obj.RefCellCornerVerts = [obj.RefCellCornerVerts;new_cells];
                obj.RefCellMidFaceVerts = [obj.RefCellMidFaceVerts;new_cells];
                obj.RefCellHigherLvls = [obj.RefCellHigherLvls;new_cells];
            end
            % Face Arrays
            if nfaces > 0
                new_faces = cell(nfaces,1);
                obj.FaceVerts = [obj.FaceVerts;new_faces];
                obj.FaceCells = [obj.FaceCells;zeros(nfaces,2)];
                obj.OrthogonalProjection = [obj.OrthogonalProjection;zeros(nfaces,2)];
                obj.FaceID = [obj.FaceID;zeros(nfaces,1)];
                obj.FaceNormal = [obj.FaceNormal;zeros(nfaces,obj.Dimension)];
                obj.FaceCenter = [obj.FaceCenter;zeros(nfaces,obj.Dimension)];
                obj.FaceArea = [obj.FaceArea;zeros(nfaces,1)];
            end
            % Vertex Arrays
            if nverts > 0
                new_verts = cell(nverts,1);
                obj.Vertices    = [obj.Vertices;zeros(nverts,obj.Dimension)]; 
                obj.VertexCells = [obj.VertexCells;new_verts];
                obj.VertexFaces = [obj.VertexFaces;new_verts];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_mesh_type(obj)
            if obj.Dimension == 2
                if length(nonzeros(obj.CellVertexNumbers)) > 1
                    obj.MeshType = 'Polygon';
                else
                    if length(obj.CellVertexNumbers) == 3
                        obj.MeshType = 'Triangle';
                    elseif length(obj.CellVertexNumbers) == 4
                        obj.MeshType = 'Quadrilateral';
                    else
                        obj.MeshType = 'Polygon';
                    end
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
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Geometry Modification Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function update_geometry_info_after_modifications(obj)
            % Update All Counts
            obj.TotalVertices = size(obj.Vertices, 1);
            obj.TotalCells = length(obj.CellVerts);
            obj.TotalFaces = length(obj.FaceVerts);
            obj.BoundaryFaces = []; obj.TotalBoundaryFaces = 0;
            obj.InteriorFaces = []; obj.TotalInteriorFaces = 0;
            % Loop through all faces
            for f=1:obj.TotalFaces
                fv = obj.FaceVerts{f}; nfv = length(fv);
                fid = obj.FaceID(f);
                if fid == 0
                    obj.TotalInteriorFaces = obj.TotalInteriorFaces + 1;
                    %obj.InteriorFaces = [obj.InteriorFaces;f];
                else
                    obj.TotalBoundaryFaces = obj.TotalBoundaryFaces + 1;
                    %obj.BoundaryFaces = [obj.BoundaryFaces;f];
                end
                fverts = obj.Vertices(fv,:);
                obj.FaceCenter(f,:) = mean(fverts);
                if obj.Dimension == 1
                    obj.FaceArea(f) = 1;
                    fc = obj.FaceCenter(f,:);
                    cc = obj.CellCenter(obj.FaceCells(f,1),:);
                    obj.FaceNormal(f,:) = (fc-cc)/abs(fc-cc);
                elseif obj.Dimension == 2
                    dxf = diff(fverts);
                    obj.FaceArea(f) = norm(dxf);
                    fc = obj.FaceCenter(f,:);
                    cc = obj.CellCenter(obj.FaceCells(f,1),:);
                    tfnorm = [dxf(2),-dxf(1)]/norm(dxf);
                    if dot(fc-cc,tfnorm) < 0, tfnorm = -1*tfnorm; end
                    obj.FaceNormal(f,:) = tfnorm;
                elseif obj.Dimension == 3
                    tfnorm = zeros(1,obj.Dimension);
                    obj.FaceArea(f) = polygonArea3d(fverts);
                    fc = obj.FaceCenter(f,:); v3 = fc;
                    for i=1:nfv
                        ii = [i,mod(i,nfv)+1];
                        v1 = fverts(ii(1),:);
                        v2 = fverts(ii(2),:);
                        dv1 = v2-v1; dv2 = v3 - v1;
                        tfnorm = tfnorm + cross(dv1,dv2)/2;
                    end
                    tfnorm = tfnorm / obj.FaceArea(f);
                    obj.FaceNormal(f,:) = tfnorm / nfv;
                end
            end
            obj.InteriorFaces = zeros(obj.TotalInteriorFaces, 1);
            obj.BoundaryFaces = zeros(obj.TotalBoundaryFaces, 1);
            % Loop through all faces again
            bfcount = 1; ifcount = 1;
            for f=1:obj.TotalFaces
              fid = obj.FaceID(f);
              if fid == 0
                obj.InteriorFaces(ifcount) = f;
                ifcount = ifcount + 1;
              else
                obj.BoundaryFaces(bfcount) = f;
                bfcount = bfcount + 1;
              end
            end
            %obj.TotalBoundaryFaces = length(obj.BoundaryFaces);
            %obj.TotalInteriorFaces = length(obj.InteriorFaces);
            % Loop through all cells
            obj.CellVertexNumbers = [];
            for c=1:obj.TotalCells
                cv = obj.CellVerts{c};
                cverts = obj.Vertices(cv,:);
                cfaces = obj.CellFaces{c};
                obj.CellCenter(c,:) = mean(cverts);
                obj.CellSurfaceArea(c) = sum(obj.FaceArea(cfaces));
                if obj.Dimension == 1
                    obj.CellVolume(c) = abs(cverts(2) - cverts(1));
                elseif obj.Dimension == 2
                    obj.CellVolume(c) = polygonArea(cverts);
                elseif obj.Dimension == 3
                    % THIS IS INCOMPLETE...
                    cc = obj.CellCenter(c,:);
                    for ff=1:length(cfaces)
                        f = cfaces(ff);
                        fc = obj.FaceCenter(f,:);
                    end
                end
            end
            obj.calculate_orthogonal_projections();
            obj.determine_cell_vertex_numbers();
            obj.determine_mesh_type();
            % Loop through all edges in 3D only
            if obj.Dimension == 3
                for e=1:obj.TotalEdges
                    ev = obj.EdgeVerts(e,:);
                    everts = obj.Vertices(ev,:);
                    obj.EdgeCenter(e,:) = mean(everts);
                    obj.EdgeLength(e) = norm(diff(everts));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calculate_orthogonal_projections(obj)
            for c=1:obj.TotalCells
                cf = obj.CellFaces{c}; ncf = length(cf); h = zeros(ncf,1);
                cv = obj.CellVerts{c}; ncv = length(cv);
                if obj.Dimension == 1
                    h(:) = obj.CellVolume(c);
                elseif obj.Dimension == 2
                    cvol = obj.CellVolume(c);
                    csa = obj.CellSurfaceArea(c);
                    for f=1:ncf
                        if ncv == 3
                            h(f) = 2*cvol/obj.FaceArea(cf(f));
                        elseif ncv == 4
                            h(f) = cvol/obj.FaceArea(cf(f));
                        else
                            if mod(ncv,2) == 0
                                h(f) = 4*cvol/csa;
                            else
                                h(f) = 2*cvol/csa + sqrt((2*cvol)/(ncv*sin(2*pi/ncv)));
                            end
                        end
                    end
                elseif obj.Dimension == 3
                    cvol = obj.CellVolume(c);
                    csa = obj.CellSurfaceArea(c);
                    for f=1:ncf
                        if ncv == 4
                            h(f) = 3*cvol/obj.FaceArea(cf(f));
                        elseif ncv == 8 && ncf == 6
                            h(f) = cvol/obj.FaceArea(cf(f));
                        else
                            h(f) = 6*cvol/csa;
                        end
                    end
                end
                for f=1:ncf
                    ff = cf(f);
                    fcells = obj.FaceCells(ff,:);
                    if fcells(1) == c
                        obj.OrthogonalProjection(ff,1) = h(f);
                    else
                        obj.OrthogonalProjection(ff,2) = h(f);
                    end
                end
%                 fcnodes = cell(ncf, 1);
%                 for f=1:ncf
%                     ff = cf(f);
%                     fverts = obj.FaceVerts{ff}; nfverts = length(fverts);
%                     tfn = zeros(1, nfverts);
%                     for i=1:nfverts
%                         for j=1:ncv
%                             if fverts(i) == cv(j);
%                                 tfn(i) = j;
%                                 break
%                             end
%                         end
%                     end
%                     fcnodes{f} = tfn;
%                 end
%                 v = obj.Vertices(cv,:);
%                 h = get_orthogonal_projection(v,fcnodes);
%                 % Place projection into array
%                 for f=1:ncf
%                     ff = cf(f);
%                     fcells = obj.FaceCells(ff,:);
%                     if fcells(1) == c
%                         obj.OrthogonalProjection(ff,1) = h(f);
%                     else
%                         obj.OrthogonalProjection(ff,2) = h(f);
%                     end
%                 end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Mesh Refinement Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function refine_mesh( obj )
            ttime = tic;
            disp('-> Begin Geometry Mesh Refinement ')
            disp(['   -> Geometry Dimension: ',num2str(obj.Dimension)])
            disp(['   -> Current Refinement Level:   ',num2str(obj.MeshRefinementLevel)])
            disp(['   -> Next Refinement Level:      ',num2str(obj.MeshRefinementLevel+1)])
            if obj.Dimension == 1
                disp(['   -> Number of Refinement Flags: ',num2str(sum(obj.CellRefinementFlag))])
                refine_1D_mesh( obj );
            elseif obj.Dimension == 2
                if strcmp(obj.OriginalMeshType, 'Triangle')
                    refine_triangle_mesh( obj );
                elseif strcmp(obj.OriginalMeshType, 'Quadrilateral')
                    refine_quad_mesh( obj );
                else
                    error('Only triangle/quad mesh types supported in 2D at this time.');
                end
            elseif obj.Dimension == 3
                error('3D mesh refinement not currently working.');
            end
            obj.MeshRefinementLevel = obj.MeshRefinementLevel + 1;
            obj.CellRefinementFlag = false(obj.TotalCells,1);
            disp(['-> Total Mesh Refinement Time:  ',num2str(toc(ttime))])
            disp(' ')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_refinement_flag(obj, c)
            obj.CellRefinementFlag(c) = true;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [next_ref_levels, ref_cells] = check_cell_refinement_differences(obj, next_ref_levels, ref_cells)
            buffer = ref_cells;
            while ~isempty(buffer)
                % get flagged cell
                c = buffer(1); buffer(1) = [];
                neighs = [];
                % Determine neighbor cells
                for i=1:length(obj.RefCellFaceCells{c})
                    neighs = [neighs,obj.RefCellFaceCells{c}{i}];
                end
                % Loop through neighbor cells and determine irregularity
                for i=1:length(neighs)
                    lvl_diff = abs(next_ref_levels(neighs(i)) - next_ref_levels(c));
                    if lvl_diff > obj.MaxIrregularity
                        next_ref_levels(neighs(i)) = next_ref_levels(neighs(i)) + 1;
                        buffer(end+1) = neighs(i);
                        ref_cells(end+1) = neighs(i);
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

