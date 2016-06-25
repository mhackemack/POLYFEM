%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          FE Handler
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
classdef FEHandler < handle
    % General FEHandler information
    properties (Access = public)
        Dimension
        Degree
        FEMType
        FEMName
        BasisType
        BasisName
        FEMLumpBool = false
    end
    % Cell-wise data structures
    properties (Access = public)
        CellQuadNodes
        CellQuadWeights
        CellBasisValues
        CellBasisGrads
        CellGradientMatrix
        CellMassMatrix
        CellStiffnessMatrix
        CellFunctionMatrix
    end
    % Face-wise data structures
    properties (Access = public)
        FaceCells
        FaceQuadNodes
        FaceQuadWeights
        FaceBasisValues
        FaceBasisGrads
        FaceMassMatrix
        FaceConformingMassMatrix
        FaceGradientMatrix
        FaceCouplingGradientMatrix
        FaceGradNormalMatrix
        FaceCouplingGradNormalMatrix
        FaceStiffnessMatrix
        FaceCouplingStiffnessMatrix
        FaceFunctionMatrix
    end
    % LD Properties
    properties (Access = public)
        LDNumCellDoF
        LDNumDoFs
        LDProjection
        LDInterpolation
        LDCellMassMatrix
        LDCellRHSMassMatrix
        LDCellStiffnessMatrix
        LDCellGradientMatrix
        LDFaceMassMatrix
        LDFaceRHSMassMatrix
        LDFaceConformingMassMatrix
        LDFaceGradientMatrix
        LDFaceRHSGradientMatrix
        LDFaceCouplingGradientMatrix
    end
    % MMS information
    properties (Access = private)
        MMSBool = 0
        GaussOrder = []
    end
    % General Mesh/DoF information
    properties (Access = private)
        MeshType
        IsOrthogonal
        CellVertexNumbers
        TotalDoFs
        TotalCells
        TotalFaces
        TotalVertices
    end
    % Function evaluation members
    properties (Access = private)
        increase_quad_degree = false
        use_only_poly_quad = false
        
        bf_cell_func
        basis_eval_func
        volume_func
        surface_func
        quadrature_func
        eval_volume_func
        eval_surface_func
        volume_bools
        surface_bools
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                               Constructor
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        function obj = FEHandler (varargin)
            n = nargin;
            global glob
            if isempty(glob), glob.print_info = true; end
            if n == 0
                % empty constructor -> do nothing
            elseif n < 3
                error('Not enough input arguments...')
            else
                if glob.print_info, disp('-> Begin FE Handler Construction.'); end
                ttime = tic; rev_str = [];
                %----------------------------------------------
                % Input Scheme:
                % 
                % 1) class     GeneralGeometry/Cartesian Geometry
                % 2) class     DoFHandler
                % 3) string    Basis Name
                % 4) logical   Lumping Boolean
                % 5) logicals  Volume Booleans 
                % 6) logicals  Surface Booleans
                %
                % optional:
                % 7) logical MMS Boolean
                % 8) integer Sub-cell gauss order
                %----------------------------------------------
                obj.Dimension     = varargin{1}.Dimension;
                obj.MeshType      = varargin{1}.MeshType;
                obj.IsOrthogonal  = varargin{1}.IsOrthogonal;
                obj.TotalCells    = varargin{1}.TotalCells;
                obj.TotalFaces    = varargin{1}.TotalFaces;
                obj.TotalVertices = varargin{1}.TotalVertices;
                obj.Degree        = varargin{2}.Degree;
                obj.FEMType       = varargin{2}.FEMType;
                obj.FEMName       = varargin{2}.FEMName;
                obj.BasisName     = varargin{3};
                obj.FEMLumpBool   = varargin{4};
                obj.volume_bools  = varargin{5};
                obj.surface_bools = varargin{6};
                if n > 6
                    obj.MMSBool = varargin{7};
                end
                if obj.MMSBool
                    if n > 7
                        if varargin{8} < 1
                            error(['Gauss quadrature less than 1, input = ',num2str(varargin{5})])
                        elseif varargin{8} >= 1
                            obj.GaussOrder = varargin{8};
                        else
                            error('Unexpected error on FEHandler Input.')
                        end
                    end
                end
                obj.set_basis_functions();
                %
                % Allocate Memory Space
                % ---------------------
                obj.allocate_memory_space();
                %
                % Set Additional LD Preliminaries if Necessary
                % --------------------------------------------
%                 if obj.BasisType == 0
%                     obj.allocate_LD_memory();
%                     for c=1:obj.TotalCells
%                         cv = varargin{2}.ConnectivityArray{c}; ncv = length(cv);
%                         verts = varargin{2}.NodeLocations(cv,:);
%                         cell_center = mean(verts); 
%                         dx = verts - ones(ncv,1)*cell_center;
%                         obj.LDProjection{c}    = zeros(obj.LDNumCellDoF, ncv);
%                         obj.LDInterpolation{c} = zeros(ncv, obj.LDNumCellDoF);
%                         % Projection Matrix
%                         obj.LDProjection{c}(1,:) = 1;
%                         obj.LDProjection{c}(2:end,:) = dx';
%                         % Interpolation Matrix
%                         obj.LDInterpolation{c}(:,1) = 1;
%                         obj.LDInterpolation{c}(:,2:end) = dx;
%                     end
%                 end
                
                
                %
                % Generate Cell-Wise Basis Space
                % --------------------------------------------------------------
                if glob.print_info, disp('   -> Begin Cell-Wise Elementary Matrix Construction.'); end
                for c=1:obj.TotalCells
                    % Print Current Cell Information
                    if glob.print_info
                        msg = sprintf('      -> Building Cell: %d of %d',c,obj.TotalCells);
                        fprintf([rev_str,msg]);
                        rev_str = repmat(sprintf('\b'), 1, length(msg));
                    end
                    % LD FEM Generation
                    % ----------------------------------------------------------
                    if varargin{2}.DoFType == 0 || obj.Degree == 0
                        cind = varargin{1}.CellVerts{c}; nvverts = length(cind);
                        verts = varargin{1}.Vertices(cind,:);
                        cfaces = varargin{1}.CellFaces{c}; nfaces = length(cfaces);
                        fcnodes = cell(nfaces,1);
                        for f=1:nfaces
                            ff = cfaces(f);
                            cfn = varargin{1}.FaceVerts{ff,1};
                            if varargin{1}.FaceCells(ff,2) == c
                                cfn = fliplr(cfn);
                            end
                            tfn = zeros(1,length(cfn));
                            for i=1:length(tfn)
                                for j=1:length(cind)
                                    if cfn(i) == cind(j);
                                        tfn(i) = j;
                                        break
                                    end
                                end
                            end
                            fcnodes{f} = tfn;
                        end
                    % All Other FEM Generation
                    % ----------------------------------------------------------
                    else
                        % Collect Face Vertex Information for Cell
                        cv = varargin{2}.CellVertexNodes{c};
                        verts = varargin{2}.NodeLocations(cv,:);
                        nvverts = length(verts);
                        cind = varargin{2}.ConnectivityArray{c};
                        cfaces = varargin{1}.CellFaces{c}; nfaces = length(cfaces);
                        fcnodes = cell(nfaces,1);
                        for f=1:nfaces
                            ff = cfaces(f);
                            if varargin{1}.FaceCells(ff,1) == c
                                cfn = varargin{2}.FaceVertexNodes{ff,1};
                            else
                                cfn = varargin{2}.FaceVertexNodes{ff,2};
                            end
                            tfn = zeros(1,length(cfn));
                            for i=1:length(tfn)
                                for j=1:length(cind)
                                    if cfn(i) == cind(j);
                                        tfn(i) = j;
                                        break
                                    end
                                end
                            end
                            fcnodes{f} = tfn;
                        end
                    end
                    % Retrieve cell-wise matrices and quadrature
%                    [MV, MS, QV, QS] = obj.bf_cell_func(nvverts,verts,fcnodes,obj.Degree,obj.FEMLumpBool,obj.volume_bools,obj.surface_bools,true,12);
                    [MV, MS, QV, QS] = obj.bf_cell_func(nvverts,verts,fcnodes,obj.Degree,obj.FEMLumpBool,obj.volume_bools,obj.surface_bools,obj.MMSBool,obj.GaussOrder);
                    % Assign volumetric matrices
                    obj.CellMassMatrix{c}      = MV{1};
                    obj.CellFunctionMatrix{c}  = MV{1}*ones(size(MV{1},1),1);
                    obj.CellStiffnessMatrix{c} = MV{2};
                    obj.CellGradientMatrix{c}  = MV{3};
                    % Assign surface matrices
                    cfnodes = varargin{2}.CellFaceNodes{c};
                    for f=1:nfaces
                        ff = cfaces(f);
                        nvf = length(cfnodes{f});
                        fcells = varargin{1}.FaceCells(ff,:);
                        if fcells(1) == c
                            obj.FaceMassMatrix{ff,1}           = MS{1}{f};
                            obj.FaceFunctionMatrix{ff,1}       = MS{1}{f}*ones(nvf, 1);
                            obj.FaceConformingMassMatrix{ff,1} = MS{1}{f}(:,varargin{2}.ConformingFaceNodeNumbering{ff,1});
                            obj.FaceGradientMatrix{ff,1}       = MS{2}{f};
%                             if obj.MMSBool || obj.Dimension == 1
%                                 obj.FaceQuadNodes{ff,1}   = QS{1}{f};
%                                 obj.FaceQuadWeights{ff,1} = QS{2}{f};
%                                 obj.FaceBasisValues{ff,1} = QS{3}{f};
%                                 obj.FaceBasisGrads{ff,1}  = QS{4}{f};
%                             end
                        elseif fcells(2) == c
                            obj.FaceMassMatrix{ff,2}           = MS{1}{f};
                            obj.FaceFunctionMatrix{ff,2}       = MS{1}{f}*ones(nvf, 1);
                            obj.FaceConformingMassMatrix{ff,2} = MS{1}{f}(:,varargin{2}.ConformingFaceNodeNumbering{ff,2});
                            obj.FaceGradientMatrix{ff,2}       = MS{2}{f};
%                             if obj.MMSBool || obj.Dimension == 1
%                                 obj.FaceQuadNodes{ff,2}   = QS{1}{f};
%                                 obj.FaceQuadWeights{ff,2} = QS{2}{f};
%                                 obj.FaceBasisValues{ff,2} = QS{3}{f};
%                                 obj.FaceBasisGrads{ff,2}  = QS{4}{f};
%                             end
                        end
                    end
                    % Assign Cell Quadrature
                    if obj.MMSBool || obj.Dimension == 1
                        obj.CellQuadNodes{c}   = QV{1};
                        obj.CellQuadWeights{c} = QV{2};
                        obj.CellBasisValues{c} = QV{3};
%                         obj.CellBasisGrads{c}  = QV{4};
                    end
                end
                % Generate Coupling Face Gradient Terms
                % -------------------------------------
                if glob.print_info, fprintf(rev_str); rev_str = []; end
                if obj.surface_bools(2) || obj.surface_bools(3) || obj.surface_bools(4)
                    if glob.print_info, disp('   -> Begin Across-Face Elementary Matrix Construction.'); end
                    for f=1:obj.TotalFaces
                        tcells = varargin{1}.FaceCells(f,:);
                        if varargin{1}.FaceID(f) == 0
                            obj.FaceCouplingGradientMatrix{f,1} = cell(obj.Dimension,1);
                            obj.FaceCouplingGradientMatrix{f,2} = cell(obj.Dimension,1);
                            cnodes1 = varargin{2}.ConnectivityArray{tcells(1)}; ncnodes1 = length(cnodes1);
                            cnodes2 = varargin{2}.ConnectivityArray{tcells(2)}; ncnodes2 = length(cnodes2);
                            fnodes1 = varargin{2}.FaceCellNodes{f,1}; %nfnodes1 = length(fnodes1);
                            fnodes2 = varargin{2}.FaceCellNodes{f,2}; %nfnodes2 = length(fnodes2);
                            [~,cind1] = fixup_node_ordering(obj.Dimension, cnodes1, fnodes1);
                            [~,cind2] = fixup_node_ordering(obj.Dimension, cnodes2, fnodes2);
                            ccind1 = varargin{2}.ConformingFaceCellNodeNumbering{f,1};
                            ccind2 = varargin{2}.ConformingFaceCellNodeNumbering{f,2};
                            for d=1:obj.Dimension
                                obj.FaceCouplingGradientMatrix{f,1}{d} = zeros(ncnodes1,ncnodes2);
                                obj.FaceCouplingGradientMatrix{f,2}{d} = zeros(ncnodes2,ncnodes1);
                                obj.FaceCouplingGradientMatrix{f,1}{d}(:,ccind2) = obj.FaceGradientMatrix{f,1}{d}(:,cind1);
                                obj.FaceCouplingGradientMatrix{f,2}{d}(:,ccind1) = obj.FaceGradientMatrix{f,2}{d}(:,cind2);
                            end
                        end
                    end
                    if varargin{1}.HasPeriodicFaces
                        for ff=1:varargin{1}.TotalBoundaryFaces
                            f = varargin{1}.BoundaryFaces(ff);
                            if varargin{1}.PeriodicBools(f)
                                tcells = varargin{1}.FaceCells(f,:);
                                oface = varargin{1}.PeriodicOppositeFaces(f);
                                ofcell = varargin{1}.PeriodicFaceCells(f);
                                obj.FaceCouplingGradientMatrix{f,1} = cell(obj.Dimension,1);
                                obj.FaceCouplingGradientMatrix{f,2} = cell(obj.Dimension,1);
                                cnodes1 = varargin{2}.ConnectivityArray{tcells(1)}; ncnodes1 = length(cnodes1);
                                cnodes2 = varargin{2}.ConnectivityArray{ofcell}; ncnodes2 = length(cnodes2);
                                fnodes1 = varargin{2}.FaceCellNodes{f,1};
                                fnodes2 = varargin{2}.FaceCellNodes{oface,1};
                                [~,cind1] = fixup_node_ordering(obj.Dimension, cnodes1, fnodes1);
                                [~,cind2] = fixup_node_ordering(obj.Dimension, cnodes2, fnodes2);
                                ccind1 = varargin{2}.ConformingFaceCellNodeNumbering{f,1};
                                ccind2 = varargin{2}.ConformingFaceCellNodeNumbering{oface,1};
                                for d=1:obj.Dimension
                                    obj.FaceCouplingGradientMatrix{f,1}{d} = zeros(ncnodes1,ncnodes2);
                                    obj.FaceCouplingGradientMatrix{f,2}{d} = zeros(ncnodes2,ncnodes1);
                                    obj.FaceCouplingGradientMatrix{f,1}{d}(:,ccind2) = obj.FaceGradientMatrix{f,1}{d}(:,cind1);
                                    obj.FaceCouplingGradientMatrix{f,2}{d}(:,ccind1) = obj.FaceGradientMatrix{oface,1}{d}(:,cind2);
                                end
                            end
                        end
                    end
                end
                obj.FaceCells = [];
                if glob.print_info
                    disp(['-> Total FE Handler Generation Time:  ',num2str(toc(ttime))])
                    disp(' ')
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                          Public Accessor Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get volume quadrature nodes by cell number
        function out = get_cell_quad_nodes(obj, cellID)
            out = obj.CellQuadNodes{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get volume quadrature weights by cell number
        function out = get_cell_quad_weights(obj, cellID)
            out = obj.CellQuadWeights{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get volume mass matrix by cell number
        function out = get_cell_mass_matrix(obj, cellID)
            out = obj.CellMassMatrix{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get volume stiffness matrix by cell number
        function out = get_cell_stiffness_matrix(obj, cellID)
            out = obj.CellStiffnessMatrix{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get volume gradient matrix by cell number
        function out = get_cell_gradient_matrix(obj, cellID)
            out = obj.CellGradientMatrix{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get volume gradient matrix by cell number by dimension number
        function out = get_cell_gradient_matrix_dim(obj, cellID, dim)
            out = obj.CellGradientMatrix{cellID}{dim};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get volume function matrix by cell number
        function out = get_cell_function_matrix(obj, cellID)
            out = obj.CellFunctionMatrix{cellID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get surface quadrature nodes by face number
        function out = get_face_quad_nodes(obj, faceID)
            out = obj.FaceQuadNodes{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get surface quadrature weights by face number
        function out = get_face_quad_weights(obj, faceID)
            out = obj.FaceQuadWeights{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get surface mass matrix by face number
        function out = get_face_mass_matrix(obj, faceID)
            out = obj.FaceMassMatrix{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get surface stiffness matrix by face number
        function out = get_face_stiffness_matrix(obj, faceID)
            out = obj.FaceStiffnessMatrix{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get surface gradient matrix by face number
        function out = get_face_gradient_matrix(obj, faceID)
            out = obj.FaceGradientMatrix{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get surface gradient matrix by face number by dimension number
        function out = get_face_gradient_matrix_dim(obj, faceID, dim)
            out = obj.FaceGradientMatrix{faceID}{dim};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get surface function matrix by face number
        function out = get_face_function_matrix(obj, faceID)
            out = obj.FaceFunctionMatrix{faceID};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get basis evaluation functor
        function out = get_basis_eval_func(obj)
            out = @obj.basis_eval_func;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                           Private Functions
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set basis function handles
        function set_basis_functions(obj)
            bbname = upper(obj.BasisName);
            if obj.Dimension == 1 && ~strcmp(bbname, 'LD')
                obj.BasisType = 1;
                obj.bf_cell_func = @bf_cell_func_1D;
                return
            elseif strcmp(bbname, 'LD')
                obj.BasisType = 0;
                obj.bf_cell_func = @bf_cell_func_LD;
                obj.basis_eval_func = @LD_basis_functions;
                return
            elseif obj.Degree == 0
                obj.BasisType = 0;
                obj.bf_cell_func = @bf_cell_func_k0;
                return
            elseif strcmp(bbname, 'LAGRANGE')
                obj.BasisType = 1;
                if strcmp(obj.MeshType,'Triangle') || strcmp(obj.MeshType,'Tetrahedron')
                    obj.bf_cell_func = @bf_cell_func_Lagrange_Simplex;
                elseif strcmp(obj.MeshType,'Quadrilateral') || strcmp(obj.MeshType,'Hexahedron')
                    obj.bf_cell_func = @bf_cell_func_Lagrange_Cartesian;
%                     if obj.IsOrthogonal
%                         obj.bf_cell_func = @bf_cell_func_Lagrange_Cartesian_Orthogonal;
%                     else
%                         obj.bf_cell_func = @bf_cell_func_Lagrange_Cartesian;
%                     end
                else
                    error('Cannot recognize mesh type.');
                end
                return
            elseif strcmp(bbname, 'SERENDIPITY')
                obj.BasisType = 2;
                if strcmp(obj.MeshType,'Triangle') || strcmp(obj.MeshType,'Tetrahedron')
                    obj.bf_cell_func = @bf_cell_func_Serendipity_Simplex;
                elseif strcmp(obj.MeshType,'Quadrilateral') || strcmp(obj.MeshType,'Hexahedron')
                    obj.bf_cell_func = @bf_cell_func_Serendipity_Cartesian;
%                     if obj.IsOrthogonal
%                         obj.bf_cell_func = @bf_cell_func_Serendipity_Cartesian_Orthogonal;
%                     else
%                         obj.bf_cell_func = @bf_cell_func_Serendipity_Cartesian;
%                     end
                else
                    error('Cannot recognize mesh type.');
                end
                return
            elseif strcmp(bbname, 'PWLD')
                obj.BasisType = 3;
                if obj.Degree == 1
                    obj.bf_cell_func = @bf_cell_func_PWLD;
                elseif obj.Degree == 2
                    obj.bf_cell_func = @bf_cell_func_PWQ;
                end
                obj.basis_eval_func = @PWLD_basis_functions;
                return
%                 obj.volume_func = @PWLD_volume;
%                 obj.surface_func = @PWLD_surface;
%                 obj.quadrature_func = @PWLD_quad_gen;
            elseif strcmp(bbname, 'PWQ')
                obj.BasisType = 3;
                obj.bf_cell_func = @bf_cell_func_PWQ;
                obj.basis_eval_func = @PWLD_basis_functions;
                return
            elseif strcmp(bbname, 'WACHSPRESS')
                obj.BasisType = 4;
                obj.bf_cell_func = @bf_cell_func_Wachspress;
                return
%                 obj.volume_func = @wachspress_volume;
%                 obj.surface_func = @wachspress_surface;
%                 obj.quadrature_func = @eval_quad_gen;
%                 obj.eval_volume_func = @wachspress_basis_functions;
%                 obj.eval_surface_func = @wachspress_basis_functions;
            elseif strcmp(bbname, 'MV')
                obj.BasisType = 4;
                obj.bf_cell_func = @bf_cell_func_MV;
                return
%                 obj.volume_func = @mean_value_volume;
%                 obj.surface_func = @mean_value_surface;
%                 obj.quadrature_func = @eval_quad_gen;
%                 obj.eval_volume_func = @mean_value_basis_functions;
%                 obj.eval_surface_func = @mean_value_basis_functions;
            elseif strcmp(bbname, 'HARMONIC')
                obj.BasisType = 4;
%                 obj.volume_func = @obj.eval_volume_ints;
%                 obj.surface_func = @obj.eval_surface_ints;
%                 obj.quadrature_func = @obj.eval_quad_gen;
            elseif strcmp(bbname, 'MAXENT') || strcmp(bbname, 'MAX_ENT')
                obj.BasisType = 4;
                obj.bf_cell_func = @bf_cell_func_max_entropy;
                obj.basis_eval_func = @max_entropy_basis_functions;
                return
%                 obj.volume_func = @max_entropy_volume;
%                 obj.surface_func = @max_entropy_surface;
%                 obj.quadrature_func = @eval_quad_gen;
%                 obj.eval_volume_func = @max_entropy_basis_functions;
%                 obj.eval_surface_func = @max_entropy_basis_functions;
            elseif strcmp(bbname, 'METRIC')
                obj.BasisType = 4;
                obj.bf_cell_func = @bf_cell_func_Metric;
                return
%                 obj.quadrature_func = @eval_quad_gen;
%                 obj.eval_volume_func = @metric_basis_functions;
%                 obj.eval_surface_func = @metric_basis_functions;
%             else
%                 slength = length(bbname);
%                 if strcmp(bbname(1:11), 'BARYCENTRIC')
%                     snum = find_dash_mark_in_string(bbname);
%                     if snum == 13
%                         obj.BasisType = 4;
%                         bename = bbname(15:end);
%                     elseif slength >= 23
%                         if strcmp(bbname(1:23), 'BARYCENTRIC SERENDIPITY')
%                             bename = bbname(27:end);
%                             if obj.Degree == 1
%                                 obj.BasisType = 4;
%                             else
%                                 obj.BasisType = 5;
%                             end
%                         else
%                             error(['Cannot read basis set name: ',obj.BasisName])
%                         end
%                     end
%                     % Assign matrix construction routines
%                     % -----------------------------------
%                     obj.volume_func = @obj.eval_volume_ints;
%                     obj.surface_func = @obj.eval_surface_ints;
%                     obj.quadrature_func = @obj.eval_quad_gen;
%                     % Assign basis function generator routines
%                     % ----------------------------------------
%                     if strcmp(bename, 'MV')
%                         obj.eval_volume_func = @mean_value_basis_functions;
%                         obj.eval_surface_func = @mean_value_basis_functions;
%                     elseif strcmp(bename, 'WACHSPRESS')
%                         obj.eval_volume_func = @wachspress_basis_functions;
%                             obj.eval_surface_func = @wachspress_basis_functions;
%                     elseif strcmp(bename, 'PWLD')
%                         obj.increase_quad_degree = true;
%                         obj.use_only_poly_quad = true;
%                         obj.eval_volume_func = @PWLD_basis_functions_volume;
%                         obj.eval_surface_func = @PWLD_basis_functions_surface;
%                     end
%                 else
%                     error(['Cannot read basis set name: ',obj.BasisName])
%                 end
            else
                error('Cannot recognize basis function type.')
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Allocate all memory space based on FEM/Basis methods
        function allocate_memory_space(obj)
            % Allocate cell-wise memory
            obj.CellMassMatrix = cell(obj.TotalCells,1);
            obj.CellStiffnessMatrix = cell(obj.TotalCells,1);
            obj.CellGradientMatrix  = cell(obj.TotalCells,1);
            obj.CellFunctionMatrix = cell(obj.TotalCells,1);
            % Allocate face-wise memory
            if obj.surface_bools(1)
                obj.FaceMassMatrix = cell(obj.TotalFaces,2);
                obj.FaceConformingMassMatrix = cell(obj.TotalFaces,2);
                obj.FaceFunctionMatrix = cell(obj.TotalFaces,2);
            end
            if obj.surface_bools(2)
                obj.FaceGradientMatrix = cell(obj.TotalFaces,2);
                obj.FaceCouplingGradientMatrix = cell(obj.TotalFaces,2);
            end
            if obj.surface_bools(3)
                
            end
            if obj.surface_bools(4)
                
            end
            obj.FaceCells = zeros(obj.TotalFaces,2);
            % Quadrature Information
            obj.CellQuadNodes = cell(obj.TotalCells,1);
            obj.CellQuadWeights = cell(obj.TotalCells,1);
            obj.CellBasisValues = cell(obj.TotalCells,1);
            obj.CellBasisGrads = cell(obj.TotalCells,1);
            obj.FaceQuadNodes = cell(obj.TotalFaces,1);
            obj.FaceQuadWeights = cell(obj.TotalFaces,1);
            obj.FaceBasisValues = cell(obj.TotalFaces,1);
            obj.FaceBasisGrads = cell(obj.TotalFaces,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Allocate LD memory space
        function allocate_LD_memory(obj)
            obj.LDNumCellDoF = obj.Dimension + 1;
            obj.LDProjection = cell(obj.TotalCells, 1);
            obj.LDInterpolation = cell(obj.TotalCells, 1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        function varargout = eval_volume_ints(obj, varargin)
            % Get Input Information
            % ---------------------
            verts = varargin{1};
            faces = varargin{2};
            flags = varargin{3};
            deg = obj.Degree + 1;
            if flags(2) || flags(3)
                grad_bool = true;
            else
                grad_bool = false;
            end
            [nv, dim] = size(verts);
            % Get Quadrature
            % --------------
            if obj.increase_quad_degree, deg = deg + 4; end
            [qx, qw] = get_general_volume_quadrature(verts, faces, deg + 2, obj.use_only_poly_quad);
            nqx = length(qw);
            % Get Basis Space
            % ---------------
            if obj.BasisType == 5
                if ~grad_bool
                    bvals = barycentric_serendipity(verts, qx, faces, obj.eval_volume_func);
                else
                    [bvals, bgrads] = barycentric_serendipity(verts, qx, faces, obj.eval_volume_func);
                end
            else
                if ~grad_bool
                    bvals = obj.eval_volume_func(verts, qx, faces);
                else
                    [bvals, bgrads] = obj.eval_volume_func(verts, qx, faces);
                end
            end
            % Set Outputs
            % -----------
            nbf = nv*obj.Degree;
            counter = 1;
            % mass matrix
            if flags(1)
                m = zeros(nbf, nbf);
                for q=1:nqx
                    bt = bvals(q,:);
                    m = m + bt'*bt*qw(q);
                end
                varargout{counter} = m;
                counter = counter + 1;
            end
            % stiffness matrix
            if flags(2)
                k = zeros(nbf, nbf);
                for q=1:nqx
                    bgt = bgrads(:,:,q);
                    k = k + bgt*bgt'*qw(q);
                end
                varargout{counter} = k;
                counter = counter + 1;
            end
            % gradient matrix
            if flags(3)
                g = cell(obj.Dimension, 1);
                for d=1:obj.Dimension
                    g{d} = zeros(nbf, nbf);
                end
                for q=1:nqx
                    bt = bvals(q,:);
                    bgt = bgrads(:,:,q);
                    for d=1:obj.Dimension
                        g{d} = g{d} + bgt(:,d)*bt*qw(q);
                    end
                end
                varargout{counter} = g;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        function varargout = eval_surface_ints(obj, varargin)
            % Get Input Information
            % ---------------------
            verts = varargin{1};
            fnodes = varargin{2};
            flags = varargin{3};
            deg = varargin{4};
            [nv, dim] = size(verts);
            % Loop through faces and compute matrices
            % ---------------------------------------
            for f=1:length(fnodes)
                % Get Quadrature
                % --------------
                [qx, qw] = get_face_quad(verts, fnodes{f}, obj.Degree+1); nqx = length(qw);
                % Get Basis Space
                % ---------------
                if flags(2)
                    grad_bool = true;
                else
                    grad_bool = false;
                end
                if obj.BasisType == 5
                    if ~grad_bool
                        bvals = barycentric_serendipity(verts, qx, fnodes{f}, obj.eval_surface_func);
                    else
                        [bvals, bgrads] = barycentric_serendipity(verts, qx, fnodes{f}, obj.eval_surface_func);
                    end
                else
                    if ~grad_bool
                        bvals = obj.eval_surface_func(verts, qx, fnodes{f});
                    else
                        [bvals, bgrads] = obj.eval_surface_func(verts, qx, fnodes{f});
                    end
                end
                % Set Outputs
                % -----------
                nds = get_surface_dofs(dim, nv, fnodes{f}, obj.Degree);
                nbf = nv*obj.Degree;
                counter = 1;
                % mass matrix
                if flags(1)
                    m = zeros(nbf, nbf);
                    for q=1:nqx
                        bt = bvals(q,:);
                        m = m + bt'*bt*qw(q);
                    end
                    varargout{counter}{f} = m(nds,nds);
                    counter = counter + 1;
                end
                % gradient matrix
                if flags(2)
                    g = cell(obj.Dimension, 1);
                    for d=1:obj.Dimension
                        g{d} = zeros(nbf, nbf);
                    end
                    for q=1:nqx
                        bt = bvals(q,:);
                        bgt = bgrads(:,:,q);
                        for d=1:obj.Dimension
                            g{d} = g{d} + bgt(:,d)*bt*qw(q);
                        end
                    end
                    varargout{counter}{f} = g;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        function varargout = eval_quad_gen(obj, varargin)
            % Get Input Information
            % ---------------------
            nout = nargout;
            verts = varargin{1};
            deg = varargin{2};
            faces = varargin{3};
            if nout > 3
                fdeg = varargin{4};
            else
                fdeg = [];
            end
            if nout > 4
                nverts = varargin{5};
            else
                nverts = [];
            end
            % Get Quadrature
            % --------------
            [nv, dim] = size(verts);
            [qx, qw] = get_general_volume_quadrature(verts, faces, deg, obj.use_only_poly_quad);
            if nout == 4
                grad_bool = true;
            else
                grad_bool = false;
            end
            if obj.BasisType == 5
                if ~grad_bool
                    bvals = barycentric_serendipity(verts, qx, faces, obj.eval_func);
                else
                    [bvals, bgrads] = barycentric_serendipity(verts, qx, faces, obj.eval_func);
                end
            else
                if ~grad_bool
                    bvals = obj.eval_func(verts, qx, faces, fdeg, nverts);
                else
                    [bvals, bgrads] = obj.eval_func(verts, qx, faces, fdeg, nverts);
                end
            end
            % Set Outputs
            % -----------
            varargout{1} = qx;
            varargout{2} = qw;
            varargout{3} = bvals;
            if grad_bool
                varargout{3} = bgrads;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Additional Function Lists - 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function out = get_surface_dofs(dim, nv, fnodes, deg)
out = fnodes;
if deg > 1
    if dim == 2
        out = [out, nv + (fnodes(1)-1)*(deg-1)+1:nv + (fnodes(1)-1)*(deg-1)+(deg-1)];
    else
        
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function out = find_dash_mark_in_string(string_in)
out = [];
for i=1:length(string_in)
    if strcmp(string_in(i), '-')
        out = i; return
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function [qx, qw] = get_face_quad(verts, fnodes, deg)
% Get Input Information
% ---------------------
[nv, dim] = size(verts);
vv = verts(fnodes,:);
% Get Quadrature Info
if dim == 1
    qx = vv; qw = 1;
elseif dim == 2
    dx = diff(verts(fnodes,:));
    flen = norm(dx);
    [qx, qw] = get_legendre_gauss_quad(deg);
    qx = [dx(1)*qx + vv(1,1), dx(2)*qx + vv(1,2)];
    qw = qw * flen;
else
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%