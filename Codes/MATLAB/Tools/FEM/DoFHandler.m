%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          DoF Handler
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef DoFHandler < handle
    properties (Access = public)
        TotalDoFs
        Degree
        FEMType         % 1 = CFEM, 2 = DGFEM, 3 = WGFEM
        FEMName
        DoFType         % 0 = LD, 1 = Lagrange, 2 = Serendipity
        MaxCellNodes
    end
    properties (Access = public)
        Dimension
        GeometryType
        
        TotalVertices
        TotalBoundaryVertices
        TotalInteriorVertices
        TotalCells
        TotalFaces
        
        VertexNodes
        BoundaryNodes
        InteriorNodes
        FaceCells
        FaceCellNodes
        CellFaceNodes
        
        ConnectivityArray
        CellVertexNodes
        FaceVertexNodes
        FaceNodePartners
        FaceCellNodeNumbering
        NodeLocations
        
        ConformingFaceNodes
        ConformingFaceNodeNumbering
        ConformingFaceCellNodeNumbering
    end
    % LD Properties
    properties (Access = public)
        LDAveragePosition
        LDNumCellDoF
        LDTotalDoFs
        LDCellDoFs
        LDCellProjection, LDCellInterpolation
        LDFaceProjection, LDFaceInterpolation
    end
    properties (Access = public)
        HasPeriodicFaces = false
        PeriodicFaceDoFs
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                           Constructor
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = DoFHandler (varargin)
            n = nargin;
            global glob
            if isempty(glob), glob.print_info = true; end
            if n == 0
                % empty constructor -> do nothing
            elseif n == 4
                % Read in Geometry and Finite Element Order
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                d = varargin{1};
                obj.Dimension = d.Dimension;
                obj.TotalVertices = d.TotalVertices;
                obj.TotalCells = d.TotalCells;
                obj.TotalFaces = d.TotalFaces;
                obj.Degree = varargin{2};
                obj.FaceCells = d.FaceCells;
                obj.VertexNodes = cell(obj.TotalVertices,1);
                obj.CellFaceNodes = cell(obj.TotalCells,1);
                obj.CellVertexNodes = cell(obj.TotalCells,1);
                obj.FaceCellNodes = cell(obj.TotalFaces,2);
                obj.FaceVertexNodes = cell(obj.TotalFaces,2);
                obj.ConnectivityArray = cell(obj.TotalCells,1);
                if strcmpi(varargin{3},'cfem')
                    obj.FEMType = 1;
                    obj.FEMName = 'CFEM';
                elseif strcmpi(varargin{3},'dfem')
                    obj.FEMType = 2;
                    obj.FEMName = 'DFEM';
                elseif strcmpi(varargin{3},'wgfem')
                    obj.FEMType = 3;
                    obj.FEMName = 'WGFEM';
                else
                    error('Unsure of FEM Type.')
                end
                obj.DoFType = varargin{4};
                if obj.DoFType==0 && obj.FEMType~=2
                    error('LD Basis Functions only supported for DGFEM.');
                end
                clear varargin
                                
                % Get Geometry Information and Build DoF Arrays
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if obj.Dimension == 1
                    obj.GeometryType = 'Simplex';
                elseif obj.Dimension == 2
                    if strcmp(d.MeshType,'Triangle')
                        obj.GeometryType = 'Simplex';
                    elseif strcmp(d.MeshType,'Quadrilateral')
                        obj.GeometryType = 'Quads';
                    else
                        obj.GeometryType = 'Polygons';
                    end
                elseif obj.Dimension == 3
                    if strcmp(d.MeshType,'Tetrahedron')
                        obj.GeometryType = 'Simplex';
                    elseif strcmp(d.MeshType,'Hexahedron')
                        obj.GeometryType = 'Quads';
                    else
                        obj.GeometryType = 'Polygons';
                    end
                else
                    error('You specified more than 3 spatial dimensions...')
                end
                ttime = tic;
                if glob.print_info, disp('-> Begin Degree of Freedom Construction.'); end
                if obj.Degree == 0
                    obj = generate0DegDoFs(obj, d);
                elseif obj.DoFType == 0
                    obj = generateLDDoFs(obj, d);
                else
                    if obj.Dimension == 1
                        obj = generate1DDoFs(obj, d);
                    elseif obj.Dimension == 2
                        obj = generate2DDoFs(obj, d);
                    elseif obj.Dimension == 3
                        obj = generate3DDoFs(obj, d);
                    end
                end
                if d.HasPeriodicFaces
                    obj = set_periodic_face_dofs( obj, d ); 
                    obj.HasPeriodicFaces = true;
                end
                obj.determine_face_node_partners();
                obj.determine_face_cell_nodes();
                obj.determine_max_cell_nodes();
                if glob.print_info
                    disp(['-> Total Degree of Freedom Generation Time:  ',num2str(toc(ttime))])
                    disp(' ')
                end
            else
                error('Wrong number of constructor arguments...')
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                         Manipulators
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_face_node_partners(obj)
            global glob
            obj.FaceNodePartners = cell(obj.TotalFaces,1);
            for f=1:obj.TotalFaces
                fcn1 = obj.FaceCellNodes{f,1};
                fcn2 = obj.FaceCellNodes{f,2};
                obj.FaceNodePartners{f} = zeros(length(fcn1),2);
                obj.FaceNodePartners{f}(:,1) = fcn1';
                if ~isempty(fcn2) % boundary
                    for i=1:length(fcn1)
                        nl1 = obj.NodeLocations(fcn1(i),:);
                        for j=1:length(fcn2)
                            nl2 = obj.NodeLocations(fcn2(j),:);
                            if norm(nl1 - nl2) < glob.small
                                obj.FaceNodePartners{f}(i,2) = fcn2(j);
                            end
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_face_cell_nodes(obj)
            obj.FaceCellNodeNumbering = cell(obj.TotalFaces, 2);
            for f=1:obj.TotalFaces
                fcells = obj.FaceCells(f,:);
                for c=1:2
                    fcn = obj.FaceCellNodes{f,c}; nfcn = length(fcn);
                    if ~isempty(fcn) && fcells(c)~=0
                        obj.FaceCellNodeNumbering{f,c} = zeros(1,nfcn);
                        fcell = fcells(c);
                        ca = obj.ConnectivityArray{fcell}; nca = length(ca);
                        for i=1:nfcn
                            for j=1:nca
                                if fcn(i) == ca(j)
                                    obj.FaceCellNodeNumbering{f,c}(i) = j;
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function determine_max_cell_nodes( obj )
            obj.MaxCellNodes = 0;
            for c=1:obj.TotalCells
                ncn = length(obj.ConnectivityArray{c});
                if ncn > obj.MaxCellNodes, obj.MaxCellNodes = ncn; end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %                            Accessors
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getTotalCells (obj)
            out = obj.TotalCells;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getDimension (obj)
            out = obj.Dimension;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getDegree (obj)
            out = obj.Degree;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getGeometryType (obj)
            out = obj.GeometryType;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getNodeLocations (obj, num)
            out = obj.NodeLocations(num,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getLocalNodes (obj, num)
            if num > obj.TotalCells
                msg = ['\nCell number: ',num2str(num),' is greater than total number of cells: ',num2str(obj.TotalCells)];
                error('ErrorTests:convertTest',msg);
            end
            out = obj.ConnectivityArray{num};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getShiftedLocalNodes (obj, stride, num)
            out = obj.ConnectivityArray{num} + stride(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getFaceCellNodes (obj, fnum, cnum)
            out = obj.FaceCellNodes{fnum,cnum};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getShiftedFaceCellNodes (obj, stride, fnum, cnum)
            out = obj.FaceCellNodes{fnum,cnum} + stride;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = getVertexNodes(obj,vert)
            out = obj.VertexNodes{vert};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generate1DDoFs(obj, mesh)
k = obj.Degree;
verts = mesh.Vertices; tnv = length(verts);
if obj.FEMType == 1
    obj.TotalDoFs = obj.TotalCells * obj.Degree + 1;
elseif obj.FEMType == 2
    obj.TotalDoFs = obj.TotalCells * (obj.Degree + 1);
end
obj.NodeLocations = zeros(obj.TotalDoFs, 1);
obj.ConformingFaceNodes = cell(mesh.TotalFaces,2);
obj.ConformingFaceNodeNumbering = cell(mesh.TotalFaces,2);
obj.ConformingFaceCellNodeNumbering = cell(mesh.TotalFaces,2);

if obj.FEMType == 1
    l = tnv;
    obj.NodeLocations(1:tnv) = verts;
    obj.BoundaryNodes = [1,2];
    obj.InteriorNodes = 3:obj.TotalDoFs;
    for c=1:obj.TotalCells
        vv = mesh.CellVerts{c};
        % higher order: k > 1
        if k > 1
            ll = (1:(k-1)) + l;
            vv = [vv, ll];
            v0 = verts(vv(1)); h = diff(verts(vv(1:2)));
            obj.NodeLocations(ll) = v0 + (1:(k-1))*(h/(k));
        end
        obj.ConnectivityArray{c} = vv;
        obj.CellVertexNodes{c} = vv(1:2);
        obj.CellFaceNodes{c}{1} = vv(1);
        obj.CellFaceNodes{c}{2} = vv(2);
        l = l + (k-1);
    end
elseif obj.FEMType == 2
    l = 0;
    for c=1:obj.TotalCells
        vv = [l+1, l+2];
        cc = mesh.CellVerts{c};
        obj.NodeLocations(vv) = verts(cc);
        % higher order: k > 1
        if k > 1
            vv = [vv, vv(end) + (1:k-1)];
            v0 = verts(cc(1)); h = diff(verts(cc));
            obj.NodeLocations(vv) = [verts(cc); (v0 + (1:(k-1))*(h/(k)))'];
        end
        obj.ConnectivityArray{c} = vv;
        obj.CellVertexNodes{c} = vv(1:2);
        obj.CellFaceNodes{c}{1} = vv(1);
        obj.CellFaceNodes{c}{2} = vv(2);
        l = l + (k+1);
    end
end
% Loop through Faces
for f=1:obj.TotalFaces
    fcells = mesh.FaceCells(f,:);
    fv = mesh.FaceVerts{f};
    fid = mesh.FaceID(f);
    % Interior Face
    if fid == 0
        cf1 =  mesh.CellFaces{fcells(1)};
        cf2 =  mesh.CellFaces{fcells(2)};
        cf11 = obj.CellFaceNodes{fcells(1)}{1};
        cf12 = obj.CellFaceNodes{fcells(1)}{2};
        cf21 = obj.CellFaceNodes{fcells(2)}{1};
        cf22 = obj.CellFaceNodes{fcells(2)}{2};
        cv1  = obj.CellVertexNodes{fcells(1)};
        cv2  = obj.CellVertexNodes{fcells(2)};
        % Cell 1
        if f == cf1(1)
            obj.FaceCellNodes{f,1} = cf11;
        elseif f == cf1(2)
            obj.FaceCellNodes{f,1} = cf12;
        end
        if obj.FaceCellNodes{f,1} == cv1(1)
            obj.ConformingFaceCellNodeNumbering{f,1} = 1;
        elseif obj.FaceCellNodes{f,1} == cv1(2)
            obj.ConformingFaceCellNodeNumbering{f,1} = 2;
        end
        % Cell 2
        if f == cf2(1)
            obj.FaceCellNodes{f,2} = cf21;
        elseif f == cf2(2)
            obj.FaceCellNodes{f,2} = cf22;
        end
        if obj.FaceCellNodes{f,2} == cv2(1)
            obj.ConformingFaceCellNodeNumbering{f,2} = 1;
        elseif obj.FaceCellNodes{f,2} == cv2(2)
            obj.ConformingFaceCellNodeNumbering{f,2} = 2;
        end
    % Boundary Face
    else
        if abs(mesh.Vertices(fv) - mesh.minX) < 1e-13
            obj.FaceCellNodes{f,1} = obj.CellFaceNodes{fcells(1)}{1};
            obj.FaceCellNodes{f,2} = obj.CellFaceNodes{fcells(1)}{1};
            obj.ConformingFaceCellNodeNumbering{f,1} = 1;
            obj.ConformingFaceCellNodeNumbering{f,2} = 1;
        elseif abs(mesh.Vertices(fv) - mesh.maxX) < 1e-13
            obj.FaceCellNodes{f,1} = obj.CellFaceNodes{fcells(1)}{2};
            obj.FaceCellNodes{f,2} = obj.CellFaceNodes{fcells(1)}{2};
            obj.ConformingFaceCellNodeNumbering{f,1} = 2;
            obj.ConformingFaceCellNodeNumbering{f,2} = 2;
        end
    end
    obj.FaceVertexNodes{f,1} = obj.FaceCellNodes{f,1};
    obj.FaceVertexNodes{f,2} = obj.FaceCellNodes{f,2};
    obj.ConformingFaceNodes{f,1} = obj.FaceCellNodes{f,1};
    obj.ConformingFaceNodes{f,2} = obj.FaceCellNodes{f,2};
    obj.ConformingFaceNodeNumbering{f,1} = 1;
    obj.ConformingFaceNodeNumbering{f,2} = 1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generate2DDoFs(obj, mesh)
global glob
verts = mesh.Vertices;    nv = size(verts,1);
cells = mesh.CellVerts;
% Determine Lagrange information if applicable
if obj.DoFType == 1
    if strcmp(mesh.MeshType,'Triangle') && obj.Degree == 3
        center_cell_bool = true;
        center_cell_num = 1;
    elseif strcmp(mesh.MeshType,'Quadrilateral') && obj.Degree > 1
        center_cell_bool = true;
        if obj.Degree == 2
            center_cell_num = 1;
        elseif obj.Degree == 3
            center_cell_num = 4;
        end
    else
        center_cell_bool = false;
        center_cell_num = 0;
    end
else
    if strcmp(mesh.MeshType,'Quadrilateral') && obj.Degree == 4
        center_cell_bool = true;
        center_cell_num = 1;
    else
        center_cell_bool = false;
        center_cell_num = 0;
    end
end
% Determine maximum degrees of freedom and allocate memory
% --------------------------------------------------------
if obj.FEMType == 1     % CFEM
    obj.TotalDoFs = nv;
    if obj.Degree > 1
        obj.TotalDoFs = obj.TotalDoFs + obj.TotalFaces*(obj.Degree - 1);
    end
    if center_cell_bool
        obj.TotalDoFs = obj.TotalDoFs + center_cell_num*mesh.TotalCells;
    end
elseif obj.FEMType == 2 % DFEM
    d = 0;
    dfem_cverts = 0;
    if obj.Degree == 1
        for i=1:obj.TotalCells
            d = d + length(cells{i});
        end
        dfem_cverts = d;
    else
        for i=1:obj.TotalCells
            nv = length(cells{i});
            dfem_cverts = dfem_cverts + nv;
            cfaces = mesh.CellFaces{i};
            d = d + nv + length(cfaces)*(obj.Degree - 1);
        end
        if center_cell_bool
            d = d + center_cell_num*mesh.TotalCells;
        end
    end
    obj.TotalDoFs = d;
end

% Generate DoF Information Sets
% -----------------------------
if glob.print_info, disp('   -> Begin Cell and Face DoF Assignment.'); end
obj.NodeLocations = zeros(obj.TotalDoFs,obj.Dimension);
obj.BoundaryNodes = [];
% CFEM
if obj.FEMType == 1
    obj.NodeLocations(1:nv,:) = verts;
    obj.ConnectivityArray = cells;
    obj.CellVertexNodes = cells;
    fcbool = ones(obj.TotalFaces,1);
    CFN = cell(obj.TotalCells,1);
    FCN = cell(obj.TotalFaces,2);
    totv = mesh.TotalVertices;
    totf = mesh.TotalFaces;
    totvf = totv+totf*(obj.Degree - 1);
    % Get Higher-Order Face Nodes
    % ---------------------------
    if obj.Degree > 1
        d = totv; nd = obj.Degree - 1;
%         obj.NodeLocations = [obj.NodeLocations;zeros(obj.TotalFaces*(obj.Degree-1), obj.Dimension)];
        for f=1:obj.TotalFaces
            fverts = mesh.Vertices(mesh.FaceVerts{f},:);
            fa = mesh.FaceArea(f);
            dx = fa / obj.Degree;
            for i=1:nd
                d = d + 1;
                obj.NodeLocations(d,:) = fverts(1,:) + (fverts(2,:) - fverts(1,:))/fa*(dx*i);
            end
        end
    end
    % Get Interior Cell DoF Information
    % ---------------------------------
    if center_cell_bool
        if center_cell_num == 1
            obj.NodeLocations(totvf+1:end,:) = mesh.CellCenter;
        else
            
        end
    end
    % Loop through cells and assign face/cell node combinations
    % ---------------------------------------------------------
    for i=1:obj.TotalCells
        f = mesh.CellFaces{i};
        CFN{i} = cell(length(f),1);
        nf = length(f);
        cvs = mesh.CellVerts{i};
        ddd = [];
        for face=1:nf
            ff = f(face);
            fvs = mesh.FaceVerts{ff};
            ffff = zeros(1,length(fvs));
            for ii=1:length(fvs)
                for jj=1:length(cvs)
                    if fvs(ii) == cvs(jj)
                        ffff(ii) = cvs(jj);
                        break
                    end
                end
            end
            obj.FaceVertexNodes{f(face),fcbool(f(face))} = ffff;
            if obj.Degree > 1
                nd = obj.Degree - 1;
                dd = totv+(ff-1)*nd+1:totv+ff*nd;
                ffff = [ffff,dd];
                ddd = [ddd,dd];
            end
            CFN{i}{face} = ffff;
            FCN{f(face),fcbool(f(face))} = ffff;
            fcbool(f(face)) = fcbool(f(face)) + 1;
        end
        if center_cell_num == 1
            ddd = [ddd,totvf+i];
        elseif center_cell_num == 4
            
        end
        obj.ConnectivityArray{i} = [obj.ConnectivityArray{i}, ddd];
    end
    obj.CellFaceNodes = CFN;
    obj.FaceCellNodes = FCN;
% DFEM
elseif obj.FEMType == 2
    % For the love of all that is holy,
    % please don't touch this...
    % Yes, I know it is not commented...
    nd = 0;
    fcbool = ones(obj.TotalFaces,1);
    fcccs = zeros(obj.TotalFaces,1);
    for i=1:obj.TotalCells
        cf = mesh.CellFaces{i}; nf = length(cf);
        cv = mesh.CellVerts{i}; nv = length(cv);
        obj.CellFaceNodes{i} = cell(nf,1);
        nodes = verts(cv,:); ncd = nv;
        obj.CellVertexNodes{i} = nd+1:nd+ncd;
        if obj.Degree > 1
            ncd = ncd + (obj.Degree - 1) * nf + center_cell_num;
        end
        if obj.Degree == 2
            nodes = [nodes;mesh.FaceCenter(cf,:)];
            ftemp_nodes = mesh.FaceCenter(cf,:);
        end
        if obj.Degree == 3
            if nv == 3
                dv12 = nodes(2,:) - nodes(1,:);
                dv23 = nodes(3,:) - nodes(2,:);
                dv31 = nodes(1,:) - nodes(3,:);
                nodes = [nodes;nodes(1,:)+dv12/3;nodes(1,:)+dv12*2/3];
                nodes = [nodes;nodes(2,:)+dv23/3;nodes(2,:)+dv23*2/3];
                nodes = [nodes;nodes(3,:)+dv31/3;nodes(3,:)+dv31*2/3];
            elseif nv == 4
                dv12 = nodes(2,:) - nodes(1,:);
                dv23 = nodes(3,:) - nodes(2,:);
                dv34 = nodes(4,:) - nodes(3,:);
                dv41 = nodes(1,:) - nodes(4,:);
                nodes = [nodes;nodes(1,:)+dv12/3;nodes(1,:)+dv12*2/3];
                nodes = [nodes;nodes(2,:)+dv23/3;nodes(2,:)+dv23*2/3];
                nodes = [nodes;nodes(3,:)+dv34/3;nodes(3,:)+dv34*2/3];
                nodes = [nodes;nodes(4,:)+dv41/3;nodes(4,:)+dv41*2/3];
            else
                for j=1:nf
                    dx = diff(nodes([j,mod(j,nf)+1],:));
                    nodes = [nodes;nodes(j,:)+dx/3;nodes(j,:)+dx*2/3];
                end
            end 
        end
        if center_cell_num == 1
            nodes = [nodes;mesh.CellCenter(i,:)];
        elseif center_cell_num == 4
            dv510 = nodes(10,:) - nodes(5,:);
            dv69  = nodes(9,:) - nodes(6,:);
            nodes = [nodes;nodes(5,:) + dv510/3];
            nodes = [nodes;nodes(6,:) + dv69/3];
            nodes = [nodes;nodes(6,:) + dv69*2/3];
            nodes = [nodes;nodes(5,:) + dv510*2/3];
        end
        for f=1:nf
            ff = cf(f);
            fv = mesh.FaceVerts{ff}; 
            nfv = length(fv);
            tfcb = fcbool(ff);
            % get indexing
            ind = zeros(1,nfv);
            cc = 1;
            for ii=1:nv
                for jj=1:nfv
                    if cv(ii) == fv(jj)
                        ind(cc) = ii;
                        cc = cc + 1;
                    end
                end
            end
            if ind(1) == 1 && ind(2) == nv
                ind = fliplr(ind);
            end
            iind = ind + nd;
            obj.FaceVertexNodes{ff,tfcb} = iind;
            if obj.Degree == 2
                ind = [ind, nv+f];
                iind = ind + nd;
            elseif obj.Degree == 3
                ind = [ind, nv+2*(f-1)+[1,2]];
                iind = ind + nd;
            end
            if fcccs(ff) == 0, fcccs(ff) = i; end
            % assign cell/face nodes
            obj.FaceCellNodes{ff,tfcb} = iind;
            obj.CellFaceNodes{i}{f} = iind;
            fcbool(ff) = fcbool(ff) + 1;
        end
        obj.ConnectivityArray{i} = nd+1:nd+ncd;
        obj.NodeLocations(nd+1:nd+ncd,:) = nodes;
        nd = nd + ncd;
    end
    for f=1:obj.TotalFaces
        fcells = mesh.FaceCells(f,:);
        if fcells(1) ~= fcccs(f) && mesh.FaceID(f) == 0
            obj.FaceCellNodes(f,:) = fliplr(obj.FaceCellNodes(f,:));
            obj.FaceVertexNodes(f,:) = fliplr(obj.FaceVertexNodes(f,:));
        end
    end
end
% Get Boundary DoFs
% -----------------
if glob.print_info, disp('   -> Begin Cell and Face DoF Cleanup Operations.'); end
for f=1:mesh.TotalBoundaryFaces
    face = mesh.BoundaryFaces(f);
    fnodes = obj.FaceCellNodes{face,1};
    obj.BoundaryNodes = [obj.BoundaryNodes,fnodes];
end
obj.BoundaryNodes = unique(obj.BoundaryNodes);
obj.InteriorNodes = 1:obj.TotalDoFs;
IB = ismember(obj.InteriorNodes,obj.BoundaryNodes);
obj.InteriorNodes(IB==1) = [];
obj.ConformingFaceNodes = cell(mesh.TotalFaces,2);
obj.ConformingFaceNodeNumbering = cell(mesh.TotalFaces,2);
obj.ConformingFaceCellNodeNumbering = cell(mesh.TotalFaces,2);
% Cleanup Face Node Numbering
for f=1:obj.TotalFaces
    fcells = mesh.FaceCells(f,:);
    fflag = mesh.FaceID(f);
    nfverts = length(mesh.FaceVerts{f});
    if fflag ~= 0, fcells(2) = []; end
    for i=1:length(fcells)
        ncverts = length(mesh.CellVerts{fcells(i)});
        fcnodes = obj.FaceCellNodes{f,i}(1:nfverts);
        cnodes = obj.ConnectivityArray{fcells(i)}(1:ncverts); ncnodes = length(cnodes);
        ff = zeros(1,length(fcnodes));
        fff = zeros(1,length(fcnodes));
        n = 1;
        % find indexing
        for ii=1:ncnodes
            for jj=1:length(fcnodes)
                if fcnodes(jj) == cnodes(ii)
                    ff(n) = cnodes(ii);
                    fff(n) = ii;
                    n = n + 1;
                    break
                end
            end
        end
        % Fix face node ordering
        if fff(2) == ncverts && fff(1) == 1
            obj.FaceCellNodes{f,i} = [ff(end:-1:1),obj.FaceCellNodes{f,i}(nfverts+1:end)];
        elseif ff(1) ~= fcnodes(1)
            obj.FaceCellNodes{f,i} = [ff,obj.FaceCellNodes{f,i}(nfverts+1:end)];
        end
%         if fff(2) == ncverts && fff(1) == 1
%             obj.FaceCellNodes{f,i} = [ff(end:-1:1),obj.FaceCellNodes{f,i}(end:-1:nfverts+1)];
%         elseif ff(1) ~= fcnodes(1)
%             obj.FaceCellNodes{f,i} = [ff,obj.FaceCellNodes{f,i}(end:-1:nfverts+1)];
%         end
        % Assign conforming face node ordering
        nnf = length(obj.FaceCellNodes{f,i});
        obj.ConformingFaceNodes{f,i} = [obj.FaceCellNodes{f,i}(nfverts:-1:1), obj.FaceCellNodes{f,i}(end:-1:nfverts+1)];
        obj.ConformingFaceNodeNumbering{f,i} = [nfverts:-1:1,nnf:-1:nfverts+1];
        cind = zeros(1,nnf); tfnodes = obj.ConformingFaceNodes{f,i};
        cnodes = obj.ConnectivityArray{fcells(i)};
        for ii=1:nnf
            for jj=1:length(cnodes)
                if tfnodes(ii) == cnodes(jj)
                    cind(ii) = jj;
                    break
                end
            end
        end
        obj.ConformingFaceCellNodeNumbering{f,i} = cind;
        % Vertex Node Numbering
        % ---------------------
        fcnodes = obj.FaceVertexNodes{f,i};
        ff = zeros(1,length(fcnodes));
        fff = zeros(1,length(fcnodes));
        n = 1;
        for ii=1:length(cnodes)
            for jj=1:length(fcnodes)
                if fcnodes(jj) == cnodes(ii)
                    ff(n) = cnodes(ii);
                    fff(n) = ii;
                    n = n + 1;
                    break
                end
            end
        end
        if fff(2) == ncverts && fff(1) == 1
            obj.FaceVertexNodes{f,i} = fliplr(ff);
        else
            obj.FaceVertexNodes{f,i} = ff;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generate3DDoFs(obj, mesh)
global glob
verts = mesh.Vertices;    nv = size(verts,1);
cells = mesh.CellVerts;
if obj.Degree > 1
    error('Cannot do higher than 1st order in 3D. I am not dealing with it.')
end
if obj.FEMType == 1     % CFEM
    obj.TotalDoFs = nv;
elseif obj.FEMType == 2 % DFEM
    d = 0;
    for i=1:obj.TotalCells
        d = d + length(cells{i});
    end
    obj.TotalDoFs = d;
end
if glob.print_info, disp('   -> Begin Cell and Face DoF Assignment.'); end
obj.NodeLocations = zeros(obj.TotalDoFs,obj.Dimension);
obj.BoundaryNodes = [];
% CFEM
if obj.FEMType == 1
    obj.NodeLocations(1:nv,:) = verts;
    obj.ConnectivityArray = cells;
    obj.CellVertexNodes = cells;
    fcbool = ones(obj.TotalFaces,1);
    CFN = cell(obj.TotalCells,1);
    for i=1:obj.TotalCells
        f = mesh.CellFaces{i};
        CFN{i} = cell(length(f),1);
        nf = length(f);
        cvs = mesh.CellVerts{i};
        ncvs = length(cvs);
        for face=1:nf
            fvs = mesh.FaceVerts{f(face)};
            ind = zeros(1,length(fvs));
            for j=1:length(fvs)
                for k=1:ncvs
                    if fvs(j) == cvs(k)
                        ind(j) = cvs(k);
                        break
                    end
                end
            end
            obj.FaceVertexNodes{f(face),fcbool(f(face))} = ind;
            CFN{i}{face} = ind;
            obj.FaceCellNodes{f(face),fcbool(f(face))} = ind;
            fcbool(f(face)) = fcbool(f(face)) + 1;
        end
    end
    obj.CellFaceNodes = CFN;
%     for i=1:obj.TotalVertices
%         obj.VertexNodes{i} = i;
%     end
% DFEM
elseif obj.FEMType == 2
    % For the love of all that is holy,
    % please don't touch this...
    % Yes, I know it is not commented...
    nd = 0;
    fcbool = ones(obj.TotalFaces,1);
    fcccs = zeros(obj.TotalFaces,1);
    for i=1:obj.TotalCells
        cf = mesh.CellFaces{i}; nf = length(cf);
        cv = mesh.CellVerts{i}; nv = length(cv);
        obj.CellFaceNodes{i} = cell(nf,1);
        nodes = verts(cv,:); ncd = nv;
        obj.CellVertexNodes{i} = nd+1:nd+ncd;
        for f=1:nf
            ff = cf(f);
            fv = mesh.FaceVerts{ff}; 
            nfv = length(fv);
            ind = zeros(1,nfv);
            for j=1:nfv
                for jj=1:nv
                    if fv(j) == cv(jj)
                        ind(j) = jj;
                        break
                    end
                end
            end
            if fcccs(ff) == 0, fcccs(ff) = i; end
            obj.FaceCellNodes{ff,fcbool(ff)} = nd*ones(1,length(ind)) + ind;
            obj.CellFaceNodes{i}{f} = nd*ones(1,length(ind)) + ind;
            fcbool(ff) = fcbool(ff) + 1;
        end
        obj.ConnectivityArray{i} = nd+1:nd+ncd;
        obj.NodeLocations(nd+1:nd+ncd,:) = nodes;
%         for j=1:length(cv)
%             obj.VertexNodes{cv(j)} = [obj.VertexNodes{cv(j)},nd+j];
%         end
        nd = nd + ncd;
    end
    for f=1:obj.TotalFaces
        fcells = mesh.FaceCells(f,:);
        if fcells(1) ~= fcccs(f) && mesh.get_face_flags(f) == 0
            obj.FaceCellNodes(f,:) = fliplr(obj.FaceCellNodes(f,:));
        end
    end
end
% Get Boundary DoFs
% -----------------
if glob.print_info, disp('   -> Begin Cell and Face DoF Cleanup Operations.'); end
for f=1:mesh.TotalBoundaryFaces
    face = mesh.BoundaryFaces(f);
    fnodes = obj.FaceCellNodes{face,1};
    obj.BoundaryNodes = [obj.BoundaryNodes,fnodes];
end
obj.BoundaryNodes = unique(obj.BoundaryNodes);
obj.InteriorNodes = 1:obj.TotalDoFs;
IB = ismember(obj.InteriorNodes,obj.BoundaryNodes);
obj.InteriorNodes(IB==1) = [];
obj.ConformingFaceNodes = cell(mesh.TotalFaces,2);
obj.ConformingFaceNodeNumbering = cell(mesh.TotalFaces,2);
% Cleanup Face Node Numbering
for f=1:obj.TotalFaces
    flag = mesh.FaceID(f);
%     fnorm = mesh.FaceNormal(f,:);
    fcells = mesh.FaceCells(f,:);
    % First Cell Ordering
    fverts = obj.FaceCellNodes{f,1}; n = length(fverts);
%     fmean = mean(obj.NodeLocations(fverts,:));
%     fnodes = [obj.NodeLocations(fverts,:);fmean];
%     plane  = createPlane(fnodes(1:3, :));
%     pnorm = [plane(4:6);plane(7:9)];
%     ppnorm = cross(pnorm(1,:), pnorm(2,:));
%     if dot(ppnorm, fnorm) > 0
%         tp = plane(7:9);
%         plane(7:9) = plane(4:6);
%         plane(4:6) = tp;
%     end
%     if dot(cross(pnorm(1,:), pnorm(2,:)), fnorm) < 0
%         plane(4:6) = -1*plane(4:6);
%     end
%     pts2d   = planePosition(fnodes, plane);
%     pt0 = pts2d(end,:); pts2d(end,:) = []; fnodes(end,:) = [];
%     pt1 = pts2d(1,:);
%     theta0  = atan2(pt1(2)-pt0(2), pt1(1)-pt0(1));
%     theta0  = mod(theta0 + 2*pi, 2*pi);
%     pts2d   = pts2d - repmat(pt0, [n 1]);
%     angle   = atan2(pts2d(:,2), pts2d(:,1));
%     angle   = mod(angle - theta0 + 4*pi, 2*pi);
%     [~, ind] = sort(angle); 
%     obj.FaceCellNodes{f,1} = fverts(ind);
    
%     fnodes = obj.NodeLocations(fverts,:);
%     xyz = fnodes - ones(n,1)*fmean;
%     w = cross(xyz(1,:),fnorm); w = w / norm(w);
%     dd = dot(w,xyz(2,:)/norm(xyz(2,:)));
%     if dd < 0
%         obj.FaceCellNodes{f,1} = fliplr(obj.FaceCellNodes{f,1});
%     end
    
    cfaces = mesh.CellFaces{fcells(1)};
    for i=1:length(cfaces)
        if f==cfaces(i)
            obj.CellFaceNodes{fcells(1)}{i} = obj.FaceCellNodes{f,1};
            break
        end
    end
    % Assign conforming face node ordering
    nnf = length(obj.FaceCellNodes{f,1});
    obj.ConformingFaceNodes{f,1} = [obj.FaceCellNodes{f,1}(n:-1:1), obj.FaceCellNodes{f,1}(end:-1:n+1)];
    obj.ConformingFaceNodeNumbering{f,1} = [n:-1:1,nnf:-1:n+1];
    cind = zeros(1,nnf); tfnodes = obj.ConformingFaceNodes{f,1};
    cnodes = obj.ConnectivityArray{fcells(1)};
    for ii=1:nnf
        for jj=1:length(cnodes)
            if tfnodes(ii) == cnodes(jj)
                cind(ii) = jj;
                break
            end
        end
    end
    obj.ConformingFaceCellNodeNumbering{f,1} = cind;
    % Second Cell Ordering if Interior Face
    if flag == 0
        fverts = obj.FaceCellNodes{f,1};
        fnodes = obj.NodeLocations(fverts,:);
        ind2 = zeros(1,n);
        fverts2 = obj.FaceCellNodes{f,2};
        fnodes2 = obj.NodeLocations(fverts2,:);
        ii=n+1;
        for i=1:n
            ii=ii-1;
            for j=1:n
                if norm(fnodes2(j,:) - fnodes(ii,:)) < 1e-14
                    ind2(i) = fverts2(j);
                    break
                end
            end
        end
        obj.FaceCellNodes{f,2} = ind2;
        
        cfaces = mesh.CellFaces{fcells(2)};
        for i=1:length(cfaces)
            if f==cfaces(i)
                obj.CellFaceNodes{fcells(2)}{i} = obj.FaceCellNodes{f,2};
                break
            end
        end
        % Assign conforming face node ordering
        nnf = length(obj.FaceCellNodes{f,2});
        obj.ConformingFaceNodes{f,2} = [obj.FaceCellNodes{f,2}(n:-1:1), obj.FaceCellNodes{f,2}(end:-1:n+1)];
        obj.ConformingFaceNodeNumbering{f,2} = [n:-1:1,nnf:-1:n+1];
        cind = zeros(1,nnf); tfnodes = obj.ConformingFaceNodes{f,2};
        cnodes = obj.ConnectivityArray{fcells(2)};
        for ii=1:nnf
            for jj=1:length(cnodes)
                if tfnodes(ii) == cnodes(jj)
                    cind(ii) = jj;
                    break
                end
            end
        end
        obj.ConformingFaceCellNodeNumbering{f,2} = cind;
    end
    % Retrieve Face Vertex Nodes
    nf = length(mesh.FaceVerts{f});
    obj.FaceVertexNodes{f,1} = obj.FaceCellNodes{f,1}(1:nf);
    if ~flag
        obj.FaceVertexNodes{f,2} = obj.FaceCellNodes{f,2}(1:nf);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = set_periodic_face_dofs( obj, mesh )
global glob
dim = mesh.Dimension;
obj.PeriodicFaceDoFs = cell(mesh.TotalFaces, 1);
pb = false(mesh.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    if mesh.PeriodicBools(ff) && ~pb(ff);
        fcenter = mesh.FaceCenter(ff,:);
        % determin face indices
        if dim == 1
            find = 1; %dbnds = [mesh.minX, mesh.maxX];
        elseif dim == 2
            if norm(fcenter(1) - mesh.minX) < glob.small || norm(fcenter(1) - mesh.maxX) < glob.small
                find = 1; %dbnds = [mesh.minX, mesh.maxX];
            elseif norm(fcenter(2) - mesh.minY) < glob.small || norm(fcenter(2) - mesh.maxY) < glob.small
                find = 2; %dbnds = [mesh.minY, mesh.maxY];
            else
                error('Could not determine boundary.');
            end
        else
            if norm(fcenter(1) - mesh.minX) < glob.small || norm(fcenter(1) - mesh.maxX) < glob.small
                find = 1; %dbnds = [mesh.minX, mesh.maxX];
            elseif norm(fcenter(2) - mesh.minY) < glob.small || norm(fcenter(2) - mesh.maxY) < glob.small
                find = 2; %dbnds = [mesh.minY, mesh.maxY];
            elseif norm(fcenter(3) - mesh.minZ) < glob.small || norm(fcenter(3) - mesh.maxZ) < glob.small
                find = 3; %dbnds = [mesh.minZ, mesh.maxZ];
            else
                error('Could not determine boundary.');
            end
        end
        ofind = 1:dim; ofind(find) = [];
        % set face information
        of = mesh.PeriodicOppositeFaces(ff);
        dfv = obj.FaceCellNodes{ff,1}; nv = length(dfv);
        ofv = obj.FaceCellNodes{of,1};
%         dfv = obj.FaceVertexNodes{ff,1}; nv = length(dfv);
%         ofv = obj.FaceVertexNodes{of,1};
        obj.PeriodicFaceDoFs{ff} = zeros(nv, 2);
        obj.PeriodicFaceDoFs{of} = zeros(nv, 2);
        for i=1:nv
            tvf = obj.NodeLocations(dfv(i), :);
            for j=1:nv
                tovf = obj.NodeLocations(ofv(j), :);
                if norm(tvf(ofind) - tovf(ofind)) < glob.small
                    obj.PeriodicFaceDoFs{ff}(i,:) = [j, ofv(j)];
                    obj.PeriodicFaceDoFs{of}(j,:) = [i, dfv(i)];
                    break
                end
            end
        end
        pb([ff,of]) = true;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generate0DegDoFs(obj, mesh)
global glob
% Retrieve some mesh information
ncells = mesh.TotalCells;
obj.TotalDoFs = ncells;
obj.NodeLocations = mesh.CellCenters;
% Loop through cells and build DoF structures
for c=1:ncells
    obj.ConnectivityArray{c} = c;
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = generateLDDoFs(obj, mesh)
global glob
ncells = mesh.TotalCells;
obj.TotalDoFs = ncells*(obj.Dimension+1);
% Loop through cells and build DoF structures
counter = 1;
for c=1:ncells
    obj.ConnectivityArray{c} = counter:(counter+obj.Dimension+1);
    
    counter = counter + (obj.Dimension+1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%