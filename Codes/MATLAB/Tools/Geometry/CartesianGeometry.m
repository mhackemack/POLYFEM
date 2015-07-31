%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Cartesian Geometry Mesh Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB class to generate all data structures necessary
%                   to fully describe a geometric domain to be used for
%                   finite element (FEM) calculations.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:          The geometry will always being as a cartesian mesh. Various
%                   routines have been created to some simple modifications to
%                   the mesh to make it unstructured. However, these
%                   unstructured mesh modification routines are mostly untested
%                   and rather simple. More functionality can easilty be added
%                   later.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef CartesianGeometry < handle
    properties (Access = public)
        Dimension
        TotalVertices
        TotalCells
        TotalFaces
        TotalEdges
        TotalInteriorFaces
        TotalBoundaryFaces
        HasPeriodicFaces = false
        IsOrthogonal = true
        IsExtruded = true
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
        x, y, z
        Lx, Ly, Lz
        minX, maxX
        minY, maxY
        minZ, maxZ
        Diameter
    end
    properties (Access = private) % AMR Variables
        OriginalCellCount
        MaxRefinementLevel
        RefinementBool
        RefinementArray
        NumberCellFlags
        CellRefinementLevel
        NextCellRefinementLevel
        CellCornerVerts
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Constructor Methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = CartesianGeometry (varargin)
            n = nargin;
            global glob
            if isempty(glob), glob.print_info = true; end
            if n == 0
                % empty constructor -> do nothing
                return
            elseif n == 1
                error('Not enough input arguments.')
            elseif n >= 2
                if glob.print_info, disp('-> Begin Cartesian Geometry Construction.'); end
                ttime = tic;
                obj.Dimension = varargin{1};
                % 1D Construction
                if obj.Dimension == 1
                    if n > 2
                        WARNING('Detected more input arguments than necessary for 1D construction.')
                        msg = ['Assuming second input argument of type "', class(varargin{2}), '" is the correct one.'];
                        WARNING(msg)
                    end
                    if isstruct(varargin{2})
                        if isfield(varargin{2}, 'x')
                            obj.x = varargin{2}.x;
                        elseif isfield(varargin{2}, 'X')
                            obj.x = varargin{2}.X;
                        elseif isfield(varargin{2}, 'verts')
                            obj.x = varargin{2}.verts;
                        elseif isfield(varargin{2}, 'Verts')
                            obj.x = varargin{2}.Verts;
                        else
                            error('Cannot determine vertex field in structure.');
                        end
                    elseif isa(varargin{2}, 'double')
                        obj.x = varargin{2};
                    end
                    % Construct 1D geometry
                    obj.Constructor_1D();
                % 2D Construction
                elseif obj.Dimension == 2
                    if n == 2
                        if ~isstruct(varargin{2})
                            error('Need struct for second input argument to hold vertex locations.')
                        else
                            xbool = 0;
                            ybool = 0;
                            if isfield(varargin{2}, 'x')
                                obj.x = varargin{2}.x;
                                xbool = 1;
                            elseif isfield(varargin{2}, 'X')
                                obj.x = varargin{2}.X;
                                xbool = 1;
                            end
                            if isfield(varargin{2}, 'y')
                                obj.y = varargin{2}.y;
                                ybool = 1;
                            elseif isfield(varargin{2}, 'Y')
                                obj.y = varargin{2}.Y;
                                ybool = 1;
                            end
                            if ~xbool || ~ybool
                                error('Insufficient fields in argument struct.')
                            end
                        end
                    elseif n == 3
                        obj.x = varargin{2};
                        obj.y = varargin{3};
                    else
                        error('Too many input arguments.')
                    end
                    % Construct 2D geometry
                    obj.Constructor_2D();
                % 3D Construction
                elseif obj.Dimension == 3
                    if n == 2
                        if ~isstruct(varargin{2})
                            error('Need struct for second input argument to hold vertex locations.')
                        else
                            xbool = 0;
                            ybool = 0;
                            zbool = 0;
                            if isfield(varargin{2}, 'x')
                                obj.x = varargin{2}.x;
                                xbool = 1;
                            elseif isfield(varargin{2}, 'X')
                                obj.x = varargin{2}.X;
                                xbool = 1;
                            end
                            if isfield(varargin{2}, 'y')
                                obj.y = varargin{2}.y;
                                ybool = 1;
                            elseif isfield(varargin{2}, 'Y')
                                obj.y = varargin{2}.Y;
                                ybool = 1;
                            end
                            if isfield(varargin{2}, 'z')
                                obj.z = varargin{2}.z;
                                zbool = 1;
                            elseif isfield(varargin{2}, 'Z')
                                obj.z = varargin{2}.Z;
                                zbool = 1;
                            end
                            if ~xbool || ~ybool || ~zbool
                                error('Insufficient fields in argument struct.')
                            end
                        end
                    elseif n == 4
                        obj.x = varargin{2};
                        obj.y = varargin{3};
                        obj.z = varargin{4};
                    else
                        error('Too many input arguments.')
                    end
                    % Construct 3D geometry
                    obj.Constructor_3D();
                end
                % Cleanup components
                % ------------------
                obj.determine_faces();
                obj.determine_vertex_components();
                obj.determine_cell_vertex_numbers();
                % Populate Initial AMR Variables
                % ------------------------------
                obj.OriginalCellCount = obj.TotalCells;
                obj.NumberCellFlags = 0;
                obj.RefinementBool = 0;
                obj.MaxRefinementLevel = 0;
                obj.RefinementArray = zeros(obj.TotalCells, 1);
                obj.CellRefinementLevel = zeros(obj.TotalCells, 1);
                obj.CellCornerVerts = ones(obj.TotalCells, 1) * (1:2^obj.Dimension);
            end
            if glob.print_info
                disp(['-> Total Cartesian Generation Time:  ',num2str(toc(ttime))])
                disp(' ')
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Methods to assist with geometry construction/manipulations
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Constructor_1D(obj)
            obj.Lx = length(obj.x);
            % Define Total Geometry Space
            obj.MeshType = 'Quads';
            obj.OriginalMeshType = 'Quads';
            if size(obj.x,1) < size(obj.x,2)
                obj.x = obj.x';
            end
            obj.Vertices = obj.x;
            obj.TotalVertices = length(obj.Vertices);
            obj.TotalCells = obj.TotalVertices - 1;
            obj.TotalFaces = obj.TotalVertices;
            obj.TotalInteriorFaces = obj.TotalFaces - 2;
            obj.TotalBoundaryFaces = 2;
            % Build Remaining Structures
            obj.Allocate_Arrays();
            for c=1:obj.TotalCells
                obj.CellVerts{c} = [c,c+1];
                obj.CellFaceVerts{c}{1} = c;
                obj.CellFaceVerts{c}{2} = c+1;
                obj.CellFaces{c} = [c,c+1];
                obj.CellVolume(c) = obj.Vertices(c+1) - obj.Vertices(c);
                obj.CellCenter(c) = (obj.Vertices(c+1) + obj.Vertices(c))/2;
                obj.FaceVerts{c} = c;
                if c==1
                    obj.CellNeighbors{c} = c+1;
                    obj.CellNeighborFaces{c} = c+1;
                    obj.CellSurfaceArea(c) = 1;
                elseif c==obj.TotalCells
                    obj.CellNeighbors{c} = c-1;
                    obj.CellNeighborFaces{c} = c;
                    obj.CellSurfaceArea(c) = 1;
                else
                    obj.CellNeighbors{c} = [c-1,c+1];
                    obj.CellNeighborFaces{c} = [c,c+1];
                    obj.CellSurfaceArea(c) = 2;
                end
            end
            obj.FaceVerts{end} = obj.TotalCells + 1;
            obj.FaceID = [1,zeros(1,obj.TotalFaces-2),1]';
            obj.FaceArea = ones(length(obj.FaceVerts),1);
            obj.FaceCenter = obj.Vertices;
            obj.FaceNormal = [-1,ones(1,obj.TotalFaces-1)]';
            for f=1:obj.TotalFaces
                if f==1
                    obj.OrthogonalProjection(1,1) = obj.x(2) - obj.x(1);
                    obj.VertexCells{1} = 1;
                    obj.FaceCells(1,1) = 1;
                elseif f==obj.TotalFaces
                    obj.OrthogonalProjection(end,1) = obj.x(end) - obj.x(end-1);
                    obj.VertexCells{obj.TotalFaces} = obj.Lx - 1;
                    obj.FaceCells(end,1) = obj.TotalCells;
                else
                    obj.OrthogonalProjection(f,1) = obj.x(f) - obj.x(f-1);
                    obj.OrthogonalProjection(f,2) = obj.x(f+1) - obj.x(f);
                    obj.VertexCells{f} = [f-1,f];
                    obj.FaceCells(f,1) = f-1;
                    obj.FaceCells(f,2) = f;
                end
            end
            obj.minX = min(obj.Vertices(:,1));
            obj.maxX = max(obj.Vertices(:,1));
            % Calculate Diameter
            obj.Diameter = obj.maxX - obj.minX;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Constructor_2D(obj)
            % Get 2D Mesh Grid
            obj.Lx = length(obj.x);
            obj.Ly = length(obj.y);
            [X,Y] = meshgrid(obj.x, obj.y);
            X = X'; Y = Y';
            % Define Total Geometry Space
            obj.MeshType = 'Quadrilateral';
            obj.OriginalMeshType = 'Quadrilateral';
            obj.Vertices = [X(:), Y(:)];
            obj.TotalVertices = size(obj.Vertices, 1);
            obj.TotalCells = (obj.Lx - 1) * (obj.Ly - 1);
            obj.TotalFaces = (obj.Lx - 1) * obj.Ly  + (obj.Ly - 1) * obj.Lx;
            obj.TotalBoundaryFaces = 2 * ((obj.Lx - 1) + (obj.Ly - 1));
            obj.TotalInteriorFaces = obj.TotalFaces - obj.TotalBoundaryFaces;
            % Build Remaining Structures
            obj.Allocate_Arrays();
            c=0; cn = (obj.Lx-1) * obj.Ly;
            for j=1:obj.Ly-1
                yy = obj.y([j,j+1]);
                for i=1:obj.Lx-1
                    c = c + 1;
                    xx = obj.x([i,i+1]);
                    obj.CellVerts{c} = [i+(j-1)*obj.Lx:i+(j-1)*obj.Lx+1,i+j*obj.Lx+1:-1:i+j*obj.Lx];
                    obj.CellFaceVerts{c}{1} = obj.CellVerts{c}([1,2]);
                    obj.CellFaceVerts{c}{2} = obj.CellVerts{c}([2,3]);
                    obj.CellFaceVerts{c}{3} = obj.CellVerts{c}([3,4]);
                    obj.CellFaceVerts{c}{4} = obj.CellVerts{c}([4,1]);
                    obj.CellVolume(c) = (xx(2) - xx(1)) * (yy(2) - yy(1));
                    obj.CellCenter(c,:) = [sum(xx)/2, sum(yy)/2];
                    obj.CellFaces{c} = [c,cn+c+(j-1)+1,c+obj.Lx-1,cn+c+(j-1)];
                    obj.CellSurfaceArea(c) = 2*abs(diff(xx)) + 2*abs(diff(yy));
                    obj.MatID(c) = 1;
                    if i==1 && j==1
                        obj.CellNeighbors{c} = [c+1,c+(obj.Lx-1)];
                        obj.CellNeighborFaces{c} = [cn+(j-1)*obj.Lx+i+1,c+(obj.Lx-1)];
                    elseif i==1 && j==(obj.Ly-1)
                        obj.CellNeighbors{c} = [c-(obj.Lx-1),c+1];
                        obj.CellNeighborFaces{c} = [c,cn+(j-1)*obj.Lx+i+1];
                    elseif i==(obj.Lx-1) && j==1
                        obj.CellNeighbors{c} = [c+(obj.Lx-1),c-1];
                        obj.CellNeighborFaces{c} = [c+(obj.Lx-1),cn+(j-1)*obj.Lx+i];
                    elseif i==(obj.Lx-1) && j==(obj.Ly-1)
                        obj.CellNeighbors{c} = [c-(obj.Lx-1),c-1];
                        obj.CellNeighborFaces{c} = [c,cn+(j-1)*obj.Lx+i];
                    elseif i==1 && j~=1 && j~=(obj.Ly-1)
                        obj.CellNeighbors{c} = [c-(obj.Lx-1),c+1,c+(obj.Lx-1)];
                        obj.CellNeighborFaces{c} = [c,cn+(j-1)*obj.Lx+i+1,c+(obj.Lx-1)];
                    elseif i==(obj.Lx-1) && j~=1 && j~=(obj.Ly-1)
                        obj.CellNeighbors{c} = [c-(obj.Lx-1),c+(obj.Lx-1),c-1];
                        obj.CellNeighborFaces{c} = [c,c+(obj.Lx-1),cn+(j-1)*obj.Lx+i];
                    elseif j==1 && i~=1 && i~=(obj.Lx-1)
                        obj.CellNeighbors{c} = [c+1,c+(obj.Lx-1),c-1];
                        obj.CellNeighborFaces{c} = [cn+(j-1)*obj.Lx+i+1,c+(obj.Lx-1),cn+(j-1)*obj.Lx+i];
                    elseif j==(obj.Ly-1) && i~=1 && i~=(obj.Lx-1)
                        obj.CellNeighbors{c} = [c-(obj.Lx-1),c+1,c-1];
                        obj.CellNeighborFaces{c} = [c,cn+(j-1)*obj.Lx+i+1,cn+(j-1)*obj.Lx+i];
                    else
                        obj.CellNeighbors{c} = [c-(obj.Lx-1),c+1,c+(obj.Lx-1),c-1];
                        obj.CellNeighborFaces{c} = [c,cn+(j-1)*obj.Lx+i+1,c+(obj.Lx-1),cn+(j-1)*obj.Lx+i];
                    end
                end
            end
            % Loop through horizontal faces
            f = 0;
            for j=1:obj.Ly
                for i=1:obj.Lx - 1
                    f = f + 1;
                    obj.FaceCenter(f,:) = [(obj.x(i) + obj.x(i+1))/2, obj.y(j)];
                    obj.FaceArea(f) = obj.x(i+1) - obj.x(i);
                    obj.FaceVerts{f} = i+(j-1)*obj.Lx:i+(j-1)*obj.Lx+1;
                    if j==1
                        obj.FaceNormal(f,:) = [0,-1];
                        obj.OrthogonalProjection(f,1) = obj.y(2) - obj.y(1);
                        obj.FaceID(f) = 1;
                        obj.FaceCells(f,1) = f;
                    elseif j==obj.Ly
                        obj.FaceNormal(f,:) = [0,1];
                        obj.OrthogonalProjection(f,1) = obj.y(end) - obj.y(end-1);
                        obj.FaceID(f) = 1;
                        obj.FaceCells(f,1) = f - (obj.Lx-1);
                        obj.FaceVerts{f} = [obj.FaceVerts{f}(2), obj.FaceVerts{f}(1)];
                    else
                        obj.FaceNormal(f,:) = [0,1];
                        obj.OrthogonalProjection(f,1) = obj.y(j) - obj.y(j-1);
                        obj.OrthogonalProjection(f,2) = obj.y(j+1) - obj.y(j);
                        obj.FaceID(f) = 0;
                        obj.FaceCells(f,:) = [f - (obj.Lx-1), f];
                        obj.FaceVerts{f} = [obj.FaceVerts{f}(2), obj.FaceVerts{f}(1)];
                    end
                end
            end
            % Loop through vertical faces
            cn = (obj.Lx-1) * obj.Ly;
            for j=1:obj.Ly - 1
                for i=1:obj.Lx
                    f = f + 1;
                    obj.FaceCenter(f,:) = [obj.x(i), (obj.y(j) + obj.y(j+1))/2];
                    obj.FaceArea(f) = obj.y(j+1) - obj.y(j);
                    obj.FaceVerts{f} = [i+(j-1)*obj.Lx,i+j*obj.Lx];
                    if i==1
                        obj.FaceNormal(f,:) = [-1,0];
                        obj.OrthogonalProjection(f,1) = obj.x(2) - obj.x(1);
                        obj.FaceID(f) = 1;
                        obj.FaceCells(f,1) = 1+(j-1)*(obj.Lx-1);
                        obj.FaceVerts{f} = [obj.FaceVerts{f}(2), obj.FaceVerts{f}(1)];
                    elseif i==obj.Lx
                        obj.FaceNormal(f,:) = [1,0];
                        obj.OrthogonalProjection(f,1) = obj.x(end) - obj.x(end-1);
                        obj.FaceID(f) = 1;
                        obj.FaceCells(f,1) = (obj.Lx-1)*j;
                    else
                        obj.FaceNormal(f,:) = [1,0];
                        obj.OrthogonalProjection(f,1) = obj.x(i) - obj.x(i-1);
                        obj.OrthogonalProjection(f,2) = obj.x(i+1) - obj.x(i);
                        obj.FaceID(f) = 0;
                        obj.FaceCells(f,:) = i+(j-1)*(obj.Lx-1)-1:i+(j-1)*(obj.Lx-1);
                    end
                end
            end
            obj.minX = min(obj.Vertices(:,1));
            obj.maxX = max(obj.Vertices(:,1));
            obj.minY = min(obj.Vertices(:,2));
            obj.maxY = max(obj.Vertices(:,2));
            % Calculate Diameter
            dx = obj.maxX - obj.minX;
            dy = obj.maxY - obj.minY;
            obj.Diameter = dx; % Default in case of equally-sized domain
            if dx >= dy
                obj.Diameter = dx;
            elseif dy >= dx
                obj.Diameter = dy;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Constructor_3D(obj)
            % Get 3D Mesh Grid
            obj.Lx = length(obj.x);
            obj.Ly = length(obj.y);
            obj.Lz = length(obj.z);
            [X, Y, Z] = meshgrid(obj.x, obj.y, obj.z);
            XX = []; YY = []; ZZ = [];
            for k=1:obj.Lz
                xt = squeeze(X(:,:,k))';
                yt = squeeze(Y(:,:,k))';
                zt = squeeze(Z(:,:,k))';
                XX = [XX; xt(:)];
                YY = [YY; yt(:)];
                ZZ = [ZZ; zt(:)];
            end
            clear X Y Z
            % Define Total Geometry Space
            obj.MeshType = 'Hexahedron';
            obj.OriginalMeshType = 'Hexahedron';
            obj.Vertices = [XX(:), YY(:), ZZ(:)];
            obj.TotalVertices = size(obj.Vertices, 1);
            obj.TotalCells = (obj.Lx - 1) * (obj.Ly - 1) * (obj.Lz - 1);
            obj.TotalFaces = (obj.Lx - 1) * (obj.Ly - 1) * obj.Lz  + ...
                             (obj.Ly - 1) * (obj.Lz - 1) * obj.Lx  + ...
                             (obj.Lz - 1) * (obj.Lx - 1) * obj.Ly;
            obj.TotalBoundaryFaces = 2 * (obj.Lx - 1) * (obj.Lz - 1) + ...
                                     2 * (obj.Ly - 1) * (obj.Lz - 1) + ...
                                     2 * (obj.Lx - 1) * (obj.Ly - 1);
            obj.TotalInteriorFaces = obj.TotalFaces - obj.TotalBoundaryFaces;
            obj.TotalEdges = ((obj.Lx-1)*obj.Ly + (obj.Ly-1)*obj.Lx)*obj.Lz + ...
                             ((obj.Lx-1)*obj.Lz + (obj.Lz-1)*obj.Lx)*obj.Ly + ...
                             ((obj.Lz-1)*obj.Ly + (obj.Ly-1)*obj.Lz)*obj.Lx;
            % Build Remaining Structures
            obj.Allocate_Arrays();
            c = 0; 
            npz = obj.Lx*obj.Ly; nnpz = (obj.Lx-1)*(obj.Ly-1);
            nlx = obj.Lz*(obj.Lx-1)*(obj.Ly-1) + (obj.Lz-1)*(obj.Lx-1)*(obj.Ly);
            for k=1:obj.Lz-1
                zz = obj.z([k,k+1]);
                for j=1:obj.Ly-1
                    yy = obj.y([j,j+1]);
                    for i=1:obj.Lx-1
                        c = c + 1;
                        xx = obj.x([i,i+1]);
                        tvec = [i+(j-1)*obj.Lx:i+(j-1)*obj.Lx+1,i+j*obj.Lx+1:-1:i+j*obj.Lx];
                        obj.CellVerts{c} = [tvec + ones(1,length(tvec))*(k-1)*npz,...
                                            tvec + ones(1,length(tvec))*(k)*npz];
                        obj.CellVolume(c) = (xx(2) - xx(1)) * (yy(2) - yy(1)) * (zz(2) - zz(1));
                        obj.CellCenter(c,:) = [sum(xx)/2, sum(yy)/2, sum(zz)/2];
                        obj.CellFaces{c} = [c,nnpz+c,nnpz*obj.Lz+i+(k-1)*obj.Ly*(obj.Lx-1)+(j-1)*(obj.Lx-1),...
                                            nnpz*obj.Lz+i+(k-1)*obj.Ly*(obj.Lx-1)+(j)*(obj.Lx-1),...
                                            nlx+(k-1)*obj.Lx*(obj.Ly-1)+(j-1)*obj.Lx+i,...
                                            nlx+(k-1)*obj.Lx*(obj.Ly-1)+(j-1)*obj.Lx+i+1];
                        obj.CellEdges{c} = [];
                        obj.MatID(c) = 1;
                    end
                end
            end
            % Loop through z-direction faces
            f = 0;
            fn = (obj.Lx-1)*(obj.Ly-1);
            for k=1:obj.Lz
                vn = (k-1)*obj.Lx*obj.Ly;
                for j=1:obj.Ly-1
                    for i=1:obj.Lx-1
                        f = f + 1;
                        obj.FaceCenter(f,:) = [(obj.x(i+1) + obj.x(i))/2, (obj.y(j+1) + obj.y(j))/2, obj.z(k)];
                        obj.FaceArea(f) = (obj.x(i+1) - obj.x(i)) * (obj.y(j+1) - obj.y(j));
                        obj.FaceVerts{f} = [vn+i+(j-1)*obj.Lx:vn+i+(j-1)*obj.Lx+1,vn+i+j*obj.Lx+1:-1:vn+i+j*obj.Lx];
                        if k==1
                            obj.OrthogonalProjection(f,1) = obj.z(2) - obj.z(1);
                            obj.FaceNormal(f,:) = [0,0,-1];
                            obj.FaceID(f) = 1;
                            obj.FaceCells(f,1) = f;
                        elseif k==obj.Lz
                            obj.OrthogonalProjection(f,1) = obj.z(end) - obj.z(end-1);
                            obj.FaceNormal(f,:) = [0,0,1];
                            obj.FaceID(f) = 1;
                            obj.FaceCells(f,1) = f-fn;
                            obj.FaceVerts{f} = fliplr(obj.FaceVerts{f});
                        else
                            obj.OrthogonalProjection(f,1) = obj.z(k) - obj.z(k-1);
                            obj.OrthogonalProjection(f,2) = obj.z(k+1) - obj.z(k);
                            obj.FaceNormal(f,:) = [0,0,1];
                            obj.FaceID(f) = 0;
                            obj.FaceCells(f,:) = [f-fn,f];
                            obj.FaceVerts{f} = fliplr(obj.FaceVerts{f});
                        end
                    end
                end
            end
            % Loop through y-direction faces
            for k=1:obj.Lz-1
                vn = (k-1)*obj.Lx*obj.Ly;
                vn2 = k*obj.Lx*obj.Ly;
                cn = (k-1)*(obj.Lx-1)*(obj.Ly-1);
                for j=1:obj.Ly
                    for i=1:obj.Lx-1
                        f = f + 1;
                        obj.FaceCenter(f,:) = [(obj.x(i+1) + obj.x(i))/2, obj.y(j), (obj.z(k+1) + obj.z(k))/2];
                        obj.FaceArea(f) = (obj.x(i+1) - obj.x(i)) * (obj.z(k+1) - obj.z(k));
                        obj.FaceVerts{f} = [vn+(j-1)*obj.Lx+i:vn+(j-1)*obj.Lx+i+1,vn2+(j-1)*obj.Lx+i+1,vn2+(j-1)*obj.Lx+i];
                        if j==1
                            obj.OrthogonalProjection(f,1) = obj.y(2) - obj.y(1);
                            obj.FaceNormal(f,:) = [0,-1,0];
                            obj.FaceID(f) = 1;
                            obj.FaceCells(f,1) = cn+i;
                            obj.FaceVerts{f} = fliplr(obj.FaceVerts{f});
                        elseif j==obj.Ly
                            obj.OrthogonalProjection(f,1) = obj.y(end) - obj.y(end-1);
                            obj.FaceNormal(f,:) = [0,1,0];
                            obj.FaceID(f) = 1;
                            obj.FaceCells(f,1) = cn+(obj.Ly-2)*(obj.Lx-1)+i;
                        else
                            obj.OrthogonalProjection(f,1) = obj.y(j) - obj.y(j-1);
                            obj.OrthogonalProjection(f,2) = obj.y(j+1) - obj.y(j);
                            obj.FaceNormal(f,:) = [0,1,0];
                            obj.FaceID(f) = 0;
                            obj.FaceCells(f,:) = [cn+(j-2)*(obj.Lx-1)+i,cn+(j-1)*(obj.Lx-1)+i];
                        end
                    end
                end
            end
            % Loop through x-direction faces
            for k=1:obj.Lz-1
                vn = (k-1)*obj.Lx*obj.Ly;
                vn2 = k*obj.Lx*obj.Ly;
                cn = (k-1)*(obj.Lx-1)*(obj.Ly-1);
                for j=1:obj.Ly-1
                    for i=1:obj.Lx
                        f = f + 1;
                        obj.FaceCenter(f,:) = [obj.x(i), (obj.y(j+1) + obj.y(j))/2, (obj.z(k+1) + obj.z(k))/2];
                        obj.FaceArea(f) = (obj.y(j+1) - obj.y(j)) * (obj.z(k+1) - obj.z(k));
                        obj.FaceVerts{f} = [vn+(j-1)*obj.Lx+i,vn+j*obj.Lx+i,vn2+j*obj.Lx+i,vn2+(j-1)*obj.Lx+i];
                        if i==1
                            obj.OrthogonalProjection(f,1) = obj.x(2) - obj.x(1);
                            obj.FaceNormal(f,:) = [-1,0,0];
                            obj.FaceID(f) = 1;
                            obj.FaceCells(f,1) = cn+(obj.Lx-1)*(j-1)+1;
                        elseif i==obj.Lx
                            obj.OrthogonalProjection(f,1) = obj.x(end) - obj.x(end-1);
                            obj.FaceNormal(f,:) = [1,0,0];
                            obj.FaceID(f) = 1;
                            obj.FaceCells(f,1) = cn+(obj.Lx-1)*j;
                            obj.FaceVerts{f} = fliplr(obj.FaceVerts{f});
                        else
                            obj.OrthogonalProjection(f,1) = obj.x(i) - obj.x(i-1);
                            obj.OrthogonalProjection(f,2) = obj.x(i+1) - obj.x(i);
                            obj.FaceNormal(f,:) = [1,0,0];
                            obj.FaceID(f) = 0;
                            obj.FaceCells(f,:) = cn+(obj.Lx-1)*(j-1)+(i-1):cn+(obj.Lx-1)*(j-1)+i;
                            obj.FaceVerts{f} = fliplr(obj.FaceVerts{f});
                        end
                    end
                end
            end
            obj.minX = min(obj.Vertices(:,1));
            obj.maxX = max(obj.Vertices(:,1));
            obj.minY = min(obj.Vertices(:,2));
            obj.maxY = max(obj.Vertices(:,2));
            obj.minZ = min(obj.Vertices(:,3));
            obj.maxZ = max(obj.Vertices(:,3));
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
            obj.FaceID = zeros(obj.TotalFaces, 1);
            obj.FaceNormal = zeros(obj.TotalFaces, obj.Dimension);
            obj.FaceCenter = zeros(obj.TotalFaces, obj.Dimension);
            obj.FaceArea = zeros(obj.TotalFaces, 1);
            % Vertex Arrays
            obj.VertexCells = cell(obj.TotalVertices, 1);
            obj.VertexFaces = cell(obj.TotalVertices, 1);
            % Edge Arrays
            if obj.Dimension == 3
                obj.EdgeVerts  = zeros(obj.TotalEdges,2);
                obj.EdgeCenter = zeros(obj.TotalEdges,3);
                obj.EdgeLength = zeros(obj.TotalEdges,1);
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
        function [] = get_2D_faces_for_edge_cals(obj)
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Geometry Modification Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function turn_2D_mesh_to_traps(obj, deg)
            % Quick Error Checks
            % ------------------
            if rem(obj.Lx-1, 2) ~= 0 || rem(obj.Ly-1, 2) ~= 0
                error('Routine requires even number of cells in x and y.')
            end
            if deg < 0 || deg >= 1
                error('Geometry modifying factor out of range: f = [0,1).')
            end
            % Loop through Levels and form Trapezoids
            % ---------------------------------------
            for i=2:(obj.Ly-1)
                if rem(i, 2) ~= 0, continue; end
                dym = obj.y(i) - obj.y(i-1);
                dyp = obj.y(i+1) - obj.y(i);
                yfact = 1;
                for j=1:obj.Lx
                    vind = (i-1)*obj.Lx + j;
                    if yfact > 0
                        dy = dyp;
                    else
                        dy = dym;
                    end
                    obj.Vertices(vind,2) = obj.Vertices(vind,2) + yfact*deg*dy;
                    yfact = -1*yfact;
                end
            end
            obj.IsOrthogonal = false;
            obj.IsExtruded = false;
            obj.update_geometry_info_after_modifications();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function randomize_all_cell_verts(obj, deg)
            % Quick Error Checks
            % ------------------
            if deg <= 0 || deg >= .5
                error('Geometry modifying factor out of range: f = (0,0.5).')
            end
            r_nums = rand(obj.TotalVertices, obj.Dimension);
            if obj.Dimension == 1
                % do nothing
            elseif obj.Dimension == 2
                for j=2:(obj.Ly-1)
                    dym = obj.y(j) - obj.y(j-1);
                    dyp = obj.y(j+1) - obj.y(j);
                    for i=2:(obj.Lx-1)
                        dxm = obj.x(i) - obj.x(i-1);
                        dxp = obj.x(i+1) - obj.x(i);
                        vind = (j-1)*obj.Lx + i;
                        xm = obj.Vertices(vind,1) - deg*dxm;
                        xp = obj.Vertices(vind,1) + deg*dxp;
                        ym = obj.Vertices(vind,2) - deg*dym;
                        yp = obj.Vertices(vind,2) + deg*dyp;
                        xx = xm + (xp-xm)*r_nums(vind,1);
                        yy = ym + (yp-ym)*r_nums(vind,2);
                        obj.Vertices(vind,:) = [xx,yy];
                    end
                end
            else
                d = 0;
                for k=2:(obj.Lz-1)
                    dzm = obj.z(k) - obj.z(k-1);
                    dzp = obj.z(k+1) - obj.z(k);
                    for j=2:(obj.Ly-1)
                        dym = obj.y(j) - obj.y(j-1);
                        dyp = obj.y(j+1) - obj.y(j);
                        for i=2:(obj.Lx-1)
                            d = d + 1;
                            dxm = obj.x(i) - obj.x(i-1);
                            dxp = obj.x(i+1) - obj.x(i);
                            vind = (k-1)*obj.Lx*obj.Ly + (j-1)*obj.Lx + i;
                            xm = obj.Vertices(vind,1) - deg*dxm;
                            xp = obj.Vertices(vind,1) + deg*dxp;
                            ym = obj.Vertices(vind,2) - deg*dym;
                            yp = obj.Vertices(vind,2) + deg*dyp;
                            zm = obj.Vertices(vind,3) - deg*dzm;
                            zp = obj.Vertices(vind,3) + deg*dzp;
                            xx = xm + (xp-xm)*r_nums(vind,1);
                            yy = ym + (yp-ym)*r_nums(vind,2);
                            zz = zm + (zp-zm)*r_nums(vind,3);
                            obj.Vertices(vind,:) = [xx,yy,zz];
                        end
                    end
                end
            end
            obj.IsOrthogonal = false;
            obj.IsExtruded = false;
            obj.update_geometry_info_after_modifications();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function update_geometry_info_after_modifications(obj)
            % Loop through all faces
            for f=1:obj.TotalFaces
                fv = obj.FaceVerts{f};
                fverts = obj.Vertices(fv,:);
                obj.FaceCenter(f,:) = mean(fverts);
                if obj.Dimension == 2
                    dxf = diff(fverts);
                    obj.FaceArea(f) = norm(dxf);
                    tfnorm = [dxf(2),-dxf(1)];
                    if norm(tfnorm - obj.FaceNormal(f,:)) < 0, tfnorm = -1*tfnorm; end
                    obj.FaceNormal(f,:) = tfnorm;
                elseif obj.Dimension == 3
                    
                end
            end
            % Loop through all cells
            for c=1:obj.TotalCells
                cverts = obj.CellVerts{c};
                cfaces = obj.CellFaces{c};
                obj.CellCenter(c,:) = mean(cverts);
                obj.CellSurfaceArea(c) = sum(obj.FaceArea(cfaces));
                if obj.Dimension == 2
                    obj.CellVolume(c) = polygonArea(cverts);
                elseif obj.Dimension == 3
                    
                end
            end
            % Loop through faces again for orthogonal projections
%             for f=1:obj.TotalFaces
%                 if obj.FaceID(f) == 0, fcells(2) = []; end
%                 for i=1:length(fcells)
%                     c = fcells(i);
%                     ncv = length(obj.CellVerts{c});
%                     ncf = length(obj.CellFaces{c});
%                     obj.OrthogonalProjection(f,i) = get_orthogonal_length(obj.Dimension, obj.CellVolume(c), ncv, obj.FaceArea(f), obj.CellSurfaceArea(c));
%                 end
%             end
            obj.calculate_orthogonal_projections();
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
        function varargout = get_mesh_grid_locs(obj)
            if ~obj.IsOrthogonal, error('Can only build mesh on orthogonal grids.'); end
            if obj.Dimension == 1
                varargout = obj.x;
            elseif obj.Dimension == 2
                if nargout ~= 2, error('2 outputs required'); end
                [xx,yy] = meshgrid(obj.x,obj.y);
                varargout{1} = xx;
                varargout{2} = yy;
            elseif obj.Dimension == 3
                if nargout ~= 3, error('3 outputs required'); end
                [xx,yy,zz] = meshgrid(obj.x,obj.y,obj.z);
                varargout{1} = xx;
                varargout{2} = yy;
                varargout{3} = zz;
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
                obj.PeriodicBools([1, obj.TotalFaces]) = true;
                obj.PeriodicOppositeFaces(1) = obj.TotalFaces;
                obj.PeriodicOppositeFaces(end) = 1;
                obj.PeriodicFaceCells(1) = obj.TotalCells;
                obj.PeriodicFaceCells(end) = 1;
                obj.PeriodicFaceVerts{1} = [1,obj.TotalVertices];
                obj.PeriodicFaceVerts{obj.TotalVertices} = [1,1];
                return
            else
                if ~isa(dim, 'char'), error('Periodic Conditions requires either x, y, or z.'); end
                dim = lower(dim);
                if strcmp( dim, 'x' )
                    ind = 1; dbnds = [obj.minX, obj.maxX];
                elseif strcmp( dim, 'y' )
                    if obj.Dimension < 1, error('No y direction in 1D.'); end
                    ind = 2; dbnds = [obj.minY, obj.maxY];
                elseif strcmp( dim, 'z' )
                    if obj.Dimension ~= 3, error('No z direction in 1D/2D.'); end
                    ind = 3; dbnds = [obj.minZ, obj.maxZ];
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
        function set_all_boundary_flags(obj, val)
            obj.FaceID(obj.BoundaryFaces) = val;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function remove_random_verices(obj, N)
            if N == 0, return; end
            rvals = randperm(obj.TotalVertices);
            vID = rvals(1:N);
            if obj.Dimension == 2
                obj.MeshType = 'Polygon';
            elseif obj.Dimension == 3
                obj.MeshType = 'Polyhedron';
            end
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
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % List of Unstructured Routines - this is not tested for AMR...
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function clear_refinement_array(obj)
%             obj.RefinementArray = zeros(obj.TotalCells, 1);
%             obj.RefinementBool = 0;
%             obj.NumberCellFlags = 0;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function set_refinement_flag(obj, cellIDs)
%             obj.RefinementBool = 1;
%             obj.RefinementArray(cellIDs) = 1;
%             obj.NumberCellFlags = obj.NumberCellFlags + length(cellIDs);
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        % function refine_mesh - performs the following actions
        % 
        % 1) 
        % 2) 
        % 3)
        % 4)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function refine_mesh(obj)
%             ttime = tic;
%             disp('-> Begin Cartesian Geometry Mesh Refinement ')
%             disp(['   -> Geometry Dimension: ',num2str(obj.Dimension)])
%             disp(['   -> Current Refinement Level:   ',num2str(obj.MaxRefinementLevel)])
%             disp(['   -> Number of Refinement Flags: ',num2str(obj.NumberCellFlags)])
%             % Get all pertinent refinement information
%             % ----------------------------------------
%             all_cells = (1:obj.TotalCells)';
%             ref_cells = find(obj.RefinementArray ~= 0);
%             obj.NextCellRefinementLevel = obj.CellRefinementLevel;
%             obj.NextCellRefinementLevel(ref_cells) = obj.NextCellRefinementLevel(ref_cells) + 1;
%             [obj.NextCellRefinementLevel, ref_cells] = obj.check_cell_refinement_differences(obj.NextCellRefinementLevel, ref_cells);
%             num_ref_cells = length(ref_cells);
%             % cells to not be refined
%             nonref_cells = setxor(all_cells, ref_cells);
%             num_nonref_cells = length(nonref_cells);
%             % further refinement info
%             next_act_ref_lvls = obj.NextCellRefinementLevel(ref_cells);
%             [sorted_next_act_ref_lvls,sort_order] = sort(next_act_ref_lvls);
%             ordered_ref_cells = ref_cells(sort_order);
%             sorted_curr_act_ref_lvls = obj.CellRefinementLevel(ordered_ref_cells);
%             new_cell_count = (obj.TotalCells - num_ref_cells) + num_ref_cells*2^obj.Dimension;
%             % Loop through flagged cells and refine
%             for i=1:length(ordered_ref_cells)
%                 c = ordered_ref_cells(i);
%                 refine_cell(obj, c);
%             end
%             % Cleanup Everything
%             % ------------------
%             obj.determine_faces();
%             obj.determine_vertex_components();
%             obj.MaxRefinementLevel = max(max(obj.CellRefinementLevel));
%             obj.clear_refinement_array();
%             if obj.Dimension == 2
%                 obj.MeshType = 'Polygon';
%             elseif obj.Dimension == 3
%                 obj.MeshType = 'Polyhedron';
%             end
%             disp(['-> Total Mesh Refinement Time:  ',num2str(toc(ttime))])
%             disp(' ')
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function [next_ref_levels, ref_cells] = check_cell_refinement_differences(obj, next_ref_levels, ref_cells)
%             buffer = ref_cells;
%             while ~isempty(buffer)
%                 % get flagged cell
%                 c = buffer(1); buffer(1) = [];
%                 neighs = obj.CellNeighbors{c};
%                 for i=1:length(neighs)
%                     lvl_diff = abs(next_ref_levels(neighs(i)) - next_ref_levels(c));
%                     switch (lvl_diff)
%                         case{0,1}
%                             % do nothing
%                         case{2}
%                             next_ref_levels(neighs(i)) = next_ref_levels(neighs(i)) + 1;
%                             buffer(end+1) = neighs(i);
%                             ref_cells(end+1) = neighs(i);
%                         otherwise
%                             % do nothing as of right now...
%                     end
%                 end
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           Accessory Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function refine_cell(obj, cellID)
% if obj.Dimension == 1
%     refine_cell_1D(obj, cellID);
% elseif obj.Dimension == 2
%     refine_cell_2D(obj, cellID);
% elseif obj.Dimension == 3
%     refine_cell_3D(obj, cellID);
% end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function refine_cell_1D(obj, cellID)
% cfaces = obj.CellFaces{cellID};
% cverts = obj.CellVertIndices{cellID};
% % Number of added mesh components
% ncells = 1;
% nverts = 1;
% nfaces = 1;
% % Number of current mesh components
% ncurf = obj.TotalFaces;
% ncurv = obj.TotalVertices;
% ncurc = obj.TotalCells;
% % Update mesh component numbers
% obj.TotalCells = ncurf + ncells;
% obj.TotalVertices = ncurv + nverts;
% obj.TotalFaces = ncurc + nfaces;
% % Add New Components
% % ------------------
% 
% end
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
    obj.OrthogonalProjection = 1e-3*obj.OrthogonalProjection;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function refine_cell_2D(obj, cellID)
% % Get Cell Refinement Information
% % -------------------------------
% ccenter = obj.CellCenter(cellID,:);
% cverts  = obj.CellVerts{cellID};
% ccverts = obj.CellCornerVerts(cellID, :);
% cfaces  = obj.CellFaces{cellID};
% cell_curr_lvl  = obj.CellRefinementLevel(cellID);
% cell_neighbors = obj.CellNeighbors{cellID};
% cnfaces = obj.CellNeighborFaces{cellID};
% cell_neighs_curr_lvl = obj.CellRefinementLevel(cell_neighbors);
% % Switch Algorithm Based on Cell Neighbor Refinement Levels
% % ---------------------------------------------------------
% diff_val_neighs_curr = int32(cell_neighs_curr_lvl) - int32(cell_curr_lvl);
% % get more cell refinement info
% num_newverts = 2^obj.Dimension - (length(cverts) - length(ccverts)) + 1;
% num_newcells = 2^obj.Dimension-1;
% newcells = obj.TotalCells+1:obj.TotalCells+num_newcells;
% % get cell boundary faces
% boundfaces = cfaces;
% bID = obj.FaceID(boundfaces); bID = (bID==0);
% boundfaces(bID) = [];
% % get cell neighbor faces
% ccn = (diff_val_neighs_curr~=0);
% cccn = (diff_val_neighs_curr==0);
% ccnfaces = cnfaces;
% ccnfaces(ccn) = [];
% ccneighs = cell_neighbors;
% ccneighs(ccn) = [];
% cccnfaces = cnfaces;
% cccnfaces(cccn) = [];
% cccneighs = cell_neighbors;
% cccneighs(cccn) = [];
% % get non-corner vertices
% non_corner_verts = 1:length(cverts);
% non_corner_verts(ismember(non_corner_verts,ccverts)) = [];
% non_corner_verts = cverts(non_corner_verts);
% num_non_corner_verts = length(non_corner_verts);
% % get faces/vertices to be split
% split_faces = [boundfaces,ccnfaces];
% num_split_faces = length(split_faces);
% split_verts = obj.FaceCenter(split_faces,:);
% [~, nvert_ord] = sort(atan2(split_verts(:,2)-ccenter(2), split_verts(:,1)-ccenter(1)));
% newverts = [split_verts(nvert_ord,:);ccenter];
% newvertind = obj.TotalVertices+1:obj.TotalVertices+num_newverts;
% split_faces = split_faces(nvert_ord);
% num_newfaces = length(split_faces) + 2^obj.Dimension;
% newfaces = obj.TotalFaces+1:obj.TotalFaces+num_newfaces;
% newintfaces = obj.TotalFaces+num_newfaces-3:obj.TotalFaces+num_newfaces;
% % update geometry info
% num_curr_faces = obj.TotalFaces;
% num_curr_verts = obj.TotalVertices;
% obj.TotalCells = obj.TotalCells + num_newcells;
% obj.TotalVertices = obj.TotalVertices + num_newverts;
% obj.TotalFaces = obj.TotalFaces + num_newfaces;
% obj.Vertices = [obj.Vertices;newverts];
% obj.allocate_more_memory(num_newverts, num_newcells, num_newfaces);
% % get some info for the new cells
% new_cell_verts = cell(2^obj.Dimension,1);
% new_cell_faces = cell(2^obj.Dimension,1);
% for c=1:2^obj.Dimension
%     ccvert_curr = cverts(ccverts(c));
%     ccvert_faces = obj.VertexFaces{ccvert_curr};
%     new_cell_verts{c} = [ccvert_curr,newvertind(end)];
%     nc_faces = ismember(cfaces,ccvert_faces);
%     new_cell_faces{c} = [new_cell_faces{c},cfaces(nc_faces)];
% end
% % Loop through faces to be split
% next_ccnfaces = [];
% for i=1:num_split_faces
%     f = split_faces(i);
%     nf = newfaces(i);
%     if obj.FaceID(f) == 0
%         next_ccnfaces = [next_ccnfaces,nf];
%     end
%     nnf = newfaces(num_split_faces+i);
%     old_face_verts = obj.FaceVerts{f};
%     for c=1:length(new_cell_faces)
%         cff = ismember(new_cell_faces{c},f); cfind = find(cff);
%         cff(~cff) = [];
%         if ~isempty(cff)
%             if cff && new_cell_verts{c}(1) == old_face_verts(1)
%                 new_cell_verts{c} = [new_cell_verts{c},newvertind(i)];
%                 new_cell_faces{c} = [new_cell_faces{c},nnf];
%             elseif cff && new_cell_verts{c}(1) == old_face_verts(2)
%                 new_cell_faces{c}(cfind) = nf;
%                 new_cell_faces{c} = [new_cell_faces{c},nnf];
%                 new_cell_verts{c} = [new_cell_verts{c},newvertind(i)];
%             end
%         end
%     end
%     obj.FaceVerts{f} = [old_face_verts(1),newvertind(i)];
%     obj.FaceArea(f) = 0.5*obj.FaceArea(f);
%     obj.FaceCenter(f,:) = mean(obj.Vertices(obj.FaceVerts{f},:));
%     % new outer face info
%     obj.FaceVerts{nf} = [newvertind(i),old_face_verts(2)];
%     obj.FaceArea(nf) = obj.FaceArea(f);
%     obj.FaceCenter(nf,:) = mean(obj.Vertices(obj.FaceVerts{nf},:));
%     obj.FaceNormal(nf,:) = obj.FaceNormal(f,:);
%     obj.FaceID(nf) = obj.FaceID(f);
%     obj.FaceCells(nf,:) = obj.FaceCells(f,:);
%     % new interior faces from splitting outer ones
%     obj.FaceVerts{nnf} = [newvertind(i),newvertind(end)];
%     obj.FaceArea(nnf) = norm(diff(obj.Vertices(obj.FaceVerts{nnf},:)));
%     obj.FaceCenter(nnf,:) = mean(obj.Vertices(obj.FaceVerts{nnf},:));
%     dx = diff(obj.Vertices(obj.FaceVerts{nnf},:)); dx=dx/norm(dx);
%     obj.FaceNormal(nnf,:)   = [-dx(2),dx(1)];
%     obj.FaceID(nnf) = 0;
%     % add some info to neighbor cells
%     tcn = cell_neighbors((ismember(cnfaces,f)));
%     if ~isempty(tcn)
%         obj.CellFaces{tcn} = [obj.CellFaces{tcn},nf];
%         obj.CellVerts{tcn} = [obj.CellVerts{tcn},newvertind(i)];
%     end
% end
% % form new interior faces from already split outer faces
% for i=1:num_non_corner_verts
%     f = num_curr_faces + 2*num_split_faces + i;
%     ncv = non_corner_verts(i);
%     obj.FaceVerts{f} = [ncv,newvertind(end)];
%     obj.FaceCenter(f,:) = mean(obj.Vertices(obj.FaceVerts{f},:));
%     obj.FaceArea(f) = norm(diff(obj.Vertices(obj.FaceVerts{f},:)));
%     dx = diff(obj.Vertices(obj.FaceVerts{f},:));
%     obj.FaceNormal(f,:)   = [-dx(2),dx(1)];
%     obj.FaceID(f) = 0;
% end
% % Loop through faces that were split
% for i=1:length(split_faces)
%     f = split_faces(i);
%     nf = newfaces(i);
%     tind = find(ismember(ccnfaces,f));
%     cn = ccneighs(tind);
%     ncnf = next_ccnfaces(tind);
%     for c=1:2^obj.Dimension
%         if c==1
%             cc = cellID;
%         else
%             cc = newcells(c-1);
%         end
%         cf = new_cell_faces{c};
%         for ff=1:length(cf)
%             if ~isempty(cn)
%                 if cf(ff) == f
%                     if obj.FaceCells(f,1) == cn
%                         obj.FaceCells(f,2) = cc;
%                     else
%                         obj.FaceCells(f,1) = cc;
%                     end
%                 elseif cf(ff) == ncnf;
%                     if obj.FaceCells(ncnf,1) == cn
%                         obj.FaceCells(ncnf,2) = cc;
%                     else
%                         obj.FaceCells(ncnf,1) = cc;
%                     end
%                 end
%             else
%                 if cf(ff) == f
%                     obj.FaceCells(f,1) = cc;
%                 elseif cf(ff) == nf
%                     obj.FaceCells(nf,1) = cc;
%                 end
%             end
%         end
%     end
% end
% % Loop through neighbor cells that were already 1 level above
% for i=1:length(cccnfaces)
%     f = cccnfaces(i);
%     fverts = obj.FaceVerts{f};
%     cn = cccneighs(i);
% %     ncn = obj.CellNeighbors{cn};
%     for c=1:2^obj.Dimension
%         if c==1
%             cc = cellID;
%         else
%             cc = newcells(c-1);
%         end
%         tvert = new_cell_verts{c}(1);
%         if fverts(1) == tvert || fverts(2) == tvert
%             if fverts(1) == tvert
%                 new_cell_verts{c} = [new_cell_verts{c},fverts(2)];
%                 new_cell_faces{c} = [new_cell_faces{c},f];
%             else
%                 new_cell_verts{c} = [new_cell_verts{c},fverts(1)];
%                 new_cell_faces{c} = [new_cell_faces{c},f];
%             end
%             if obj.FaceCells(f,1) == cn
%                 obj.FaceCells(f,2) = cc;
%             else
%                 obj.FaceCells(f,1) = cc;
%             end
%         end
%     end
% end
% % Loop through the new interior faces
% for i=1:length(newintfaces)
%     f = newintfaces(i);
%     for c=1:2^obj.Dimension
%         if c==1
%             cc = cellID;
%         else
%             cc = newcells(c-1);
%         end
%         fverts = obj.FaceVerts{f};
%         tcverts = new_cell_verts{c};
% %         cvind = ismember(cverts,fverts);
% %         cvind(~cvind) = [];
%         cvind = logical(zeros(length(fverts),1));
%         for j=1:length(fverts)
%             for k=1:length(tcverts)
%                 if fverts(j) == tcverts(k)
%                     cvind(j) = true;
%                     break
%                 end
%             end
%         end
%         if sum(cvind) == length(fverts)
%             if obj.FaceCells(f,1) == 0
%                 obj.FaceCells(f,1) = cc;
%             else
%                 obj.FaceCells(f,2) = cc;
%             end
%             new_cell_faces{c} = [new_cell_faces{c},f];
%         end
%     end
% end
% % remove duplicate faces
% for c=1:2^obj.Dimension
%     new_cell_faces{c} = unique(new_cell_faces{c});
% end
% % loop through cells and fixup 
% for c=1:2^obj.Dimension
%     if c==1
%         cc = cellID;
%     else
%         cc = newcells(c-1);
%     end
%     obj.CellVerts{cc} = new_cell_verts{c};
%     tcverts = new_cell_verts{c};
%     ccverts = obj.Vertices(new_cell_verts{c},:);
%     cvmean = mean(ccverts);
%     [~,vord] = sort(atan2(ccverts(:,2)-cvmean(2), ccverts(:,1)-cvmean(1)));
%     obj.CellVerts{cc} = tcverts(vord);
%     obj.CellCornerVerts(cc,:) = 1:2^obj.Dimension;
%     obj.CellCenter(cc,:) = cvmean;
%     obj.CellVolume(cc) = polygonArea(obj.Vertices(obj.CellVerts{cc},:));
%     obj.CellFaces{cc} = new_cell_faces{c};
%     obj.CellSurfaceArea(cc) = sum(obj.FaceArea(obj.CellFaces{cc}));
%     % update neighbors
%     obj.CellNeighbors{cc} = [];
%     obj.CellNeighborFaces{cc} = [];
%     for i=1:length(obj.CellFaces{cc})
%         f = obj.CellFaces{cc}(i);
%         if obj.FaceID(f) == 0
%             fc = obj.FaceCells(f,:);
%             if fc(1) == cc
%                 obj.CellNeighbors{cc} = [obj.CellNeighbors{cc},fc(2)];
%             elseif fc(2) == cc
%                 obj.CellNeighbors{cc} = [obj.CellNeighbors{cc},fc(1)];
%             end
%             obj.CellNeighborFaces{cc} = [obj.CellNeighborFaces{cc},f];
%         end
%     end
% end
% % Loop through neighbor cells that gained a vertex
% for i=1:length(ccneighs)
%     c = ccneighs(i);
%     tcverts = obj.CellVerts{c};
%     ccold = obj.CellVerts{c}(obj.CellCornerVerts(c,:));
%     ccverts = obj.Vertices(tcverts,:);
%     cvmean = mean(ccverts);
%     [~,vord] = sort(atan2(ccverts(:,2)-cvmean(2), ccverts(:,1)-cvmean(1)));
%     obj.CellVerts{c} = tcverts(vord);
%     obj.CellCornerVerts(c,:) = find(ismember(tcverts(vord),ccold));
%     obj.CellNeighbors{c} = [];
%     obj.CellNeighborFaces{c} = [];
%     cf = obj.CellFaces{c};
%     for ii=1:length(cf)
%         f = cf(ii);
%         if obj.FaceID(f) == 0
%             fc = obj.FaceCells(f,:);
%             if fc(1) == c
%                 obj.CellNeighbors{c} = [obj.CellNeighbors{c},fc(2)];
%             else
%                 obj.CellNeighbors{c} = [obj.CellNeighbors{c},fc(1)];
%             end
%             obj.CellNeighborFaces{c} = [obj.CellNeighborFaces{c},f];
%         end
%     end
% end
% % update orthogonal projections and cell normals
% loop_faces = [cfaces,newfaces];
% for i=1:length(loop_faces)
%     f = loop_faces(i);
%     fcells = obj.FaceCells(f,:);
%     if obj.FaceID(f) ~= 0
%         fcells(2) = [];
%     else
%         ccenter = obj.CellCenter(fcells(1),:);
%         if dot(obj.FaceNormal(f,:), obj.FaceCenter(f,:) - ccenter) < 0
%             obj.FaceCells(f,:) = fliplr(obj.FaceCells(f,:));
%         end
%     end
%     for c=1:length(fcells)
%         cc = fcells(c);
%         cv = obj.CellVerts{cc};
%         obj.OrthogonalProjection(f,c) = get_orthogonal_length(obj.Dimension,obj.CellVolume(cc),length(cv),obj.FaceArea(f),obj.CellSurfaceArea(cc));
%     end
% end
% % Increment cell refinement level
% obj.CellRefinementLevel(cellID) = obj.CellRefinementLevel(cellID) + 1;
% obj.CellRefinementLevel(newcells) = obj.CellRefinementLevel(cellID);
% obj.MatID(newcells) = obj.MatID(cellID);
% % obj.determine_vertex_components();
% % Update Vertex Components
% verts_to_check = [cverts,newvertind];
% cells_to_check = [cellID,cell_neighbors,newcells];
% faces_to_check = [cfaces,newfaces];
% for i=1:length(cell_neighbors)
%     faces_to_check = [faces_to_check,obj.CellFaces{cell_neighbors(i)}];
% end
% faces_to_check = unique(faces_to_check);
% for i=1:length(verts_to_check)
%     obj.VertexCells{verts_to_check(i)} = [];
%     obj.VertexFaces{verts_to_check(i)} = [];
% end
% for c=1:length(cells_to_check)
%     cc = cells_to_check(c);
%     cv = obj.CellVerts{cc};
%     for i=1:length(cv)
%         for j=1:length(verts_to_check)
%             if cv(i) == verts_to_check(j)
%                 obj.VertexCells{verts_to_check(j)} = [obj.VertexCells{verts_to_check(j)},cc];
%             end
%         end
%     end
% end
% for f=1:length(faces_to_check)
%     ff = faces_to_check(f);
%     fv = obj.FaceVerts{ff};
%     for i=1:length(fv)
%         for j=1:length(verts_to_check)
%             if fv(i) == verts_to_check(j)
%                 obj.VertexFaces{verts_to_check(j)} = [obj.VertexFaces{verts_to_check(j)},ff];
%             end
%         end
%     end
% end
% for i=1:length(verts_to_check)
%     obj.VertexCells{verts_to_check(i)} = unique(obj.VertexCells{verts_to_check(i)});
%     obj.VertexFaces{verts_to_check(i)} = unique(obj.VertexFaces{verts_to_check(i)});
% end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function refine_cell_3D(obj, cellID)
% % Get Cell Refinement Information
% % -------------------------------
% cverts = obj.CellVerts{cellID};
% ccverts = obj.CellCornerVerts(cellID, :);
% cell_curr_lvl = obj.CellRefinementLevel(cellID);
% cell_next_lvl = obj.NextCellRefinementLevel(cellID);
% cell_neighbors = obj.CellNeighbors{cellID};
% cell_neighs_curr_lvl = obj.CellRefinementLevel(cell_neighbors);
% cell_neighs_next_lvl = obj.NextCellRefinementLevel(cell_neighbors);
% % Switch Algorithm Based on Cell Neighbor Refinement Levels
% % ---------------------------------------------------------
% diff_val_neighs_curr = int32(cell_neighs_curr_lvl) - int32(cell_curr_lvl);
% diff_val_neighs_next = int32(cell_neighs_next_lvl) - int32(cell_next_lvl);
% % Loop through cell neighbors
% for i = 1:length(cell_neighbors)
%     cn = cell_neighbors(i);
%     cnface = obj.CellNeighborFaces{cellID}(i);
%     cnverts = obj.CellVerts{cn};
%     cncverts = obj.CellCornerVerts(cn, :);
%     diff_curr = diff_val_neighs_curr(i);
%     diff_next = diff_val_neighs_next(i);
%     if diff_curr == 0
%         
%     elseif diff_curr == 1
%         
%     elseif diff_curr == -1
%         error('Refinement ordering gone bad...')
%     else
%         error(['Too many level differences between cell ',num2str(cellID),' and neighbor ',num2str(cn)])
%     end
% end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function refine_cell_2D_rev1(obj, cellID)
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%