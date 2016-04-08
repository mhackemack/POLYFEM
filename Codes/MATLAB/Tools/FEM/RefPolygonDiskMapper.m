%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Polygon-to-disk Mapper
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:         1) Number of Polygon Vertices
%                   2) Basis Function Name
%                   3) Basis Function Order
%                   4) Quadrature Order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef RefPolygonDiskMapper < handle
    % General Mapper Variables
    properties (Access = public)
        NumberVertices
        NumberBasisFunctions
        BasisFunctionName
        BasisFunctionOrder
        QuadratureOrder
    end
    % Reference Polygon Variables
    properties (Access = public)
        PolyVertices
        PolyFaceNodes
        ImagPolyNodes
        ImagPoly
        CellQuadNodes
        CellQuadWeights
        CellBasisValues
        CellBasisGrads
        FaceQuadNodes
        FaceQuadWeights
        FaceBasisValues
        FaceBasisGrads
    end
    % Unit Disk Variables
    properties (Access = public)
        DiskMap
        
    end
    % Reference Triangle Variables
    properties (Access = private)
        NumberRefTriNodes
        RefTriNodes
        RefTriWeights
    end
    % Unit Line Variables
    properties (Access = private)
        NumberLineNodes
        RefLineNodes
        RefLineWeights
    end
    % Basis Variables
    properties (Access = private)
        BasisFunc
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                               Constructor
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        function obj = RefPolygonDiskMapper (varargin)
            % Read in input arguments
            % ------------------------------------------------------------------
            obj.NumberVertices       = varargin{1};
            obj.BasisFunctionName    = varargin{2};
            obj.BasisFunctionOrder   = varargin{3};
            obj.QuadratureOrder      = varargin{4};
            obj.NumberBasisFunctions = obj.NumberVertices*obj.BasisFunctionOrder;
            % Build Unit Line
            % ------------------------------------------------------------------
            [obj.RefLineNodes, obj.RefLineWeights] = get_legendre_gauss_quad(obj.QuadratureOrder);
            obj.NumberLineNodes = length(obj.RefLineWeights);
            % Build Reference Triangle
            % ------------------------------------------------------------------
            [obj.RefTriNodes, obj.RefTriWeights] = TriGaussPoints(obj.QuadratureOrder);
            obj.NumberRefTriNodes = length(obj.RefTriWeights);
            obj.RefTriWeights = obj.RefTriWeights / sum(obj.RefTriWeights);
            % Build Reference Polygon
            % ------------------------------------------------------------------
            [obj.PolyVertices,obj.PolyFaceNodes] = RegularPolygon(obj.NumberVertices,1/2);
            obj.PolyFaceNodes = cell(obj.NumberVertices,1);
            for f=1:obj.NumberVertices
                obj.PolyFaceNodes{f} = [f,mod(f,obj.NumberVertices)+1];
            end
            % Convert to Imaginary Polygon and Generate SCCM Disk
            % ------------------------------------------------------------------
            obj.ImagPolyNodes = obj.PolyVertices(:,1) + 1i*obj.PolyVertices(:,2);
            obj.ImagPoly = polygon(obj.ImagPolyNodes);
            % SCCM build options
            opt.TraceSolution = 'off'; opt.Tolerance = 1e-10;
            opt.SolverMethod = 'trust'; opt.InitialGuess = [];
            % Build disk mapper
            obj.DiskMap = diskmap(obj.ImagPoly, opt);
            obj.DiskMap = center(obj.DiskMap, 0+0i);
            % Build Reference Polygon Quadrature and Basis Functions
            % ------------------------------------------------------------------
            obj.generate_volume_quadrature();
            obj.generate_surface_quadrature();
            obj.calculate_refpoly_bfgs();
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                          Public Accessor Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                           Private Functions
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set basis function generation functor
        function set_basis_functions(obj)
            if strcmpi(obj.BasisFunctionName, 'LD')
                obj.BasisFunc = @LD_basis_functions;
            elseif strcmpi(obj.BasisFunctionName, 'Wachspress')
                obj.BasisFunc = @wachspress_basis_functions;
            elseif strcmpi(obj.BasisFunctionName, 'MV')
                obj.BasisFunc = @mean_value_basis_functions;
            elseif strcmpi(obj.BasisFunctionName, 'MAXENT')
                obj.BasisFunc = @max_entropy_basis_functions;
            elseif strcmpi(obj.BasisFunctionName, 'METRIC')
                obj.BasisFunc = @metric_basis_functions;
            elseif strcmpi(obj.BasisFunctionName, 'PWLD')
                obj.BasisFunc = @PWLD_basis_functions;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Build all volume quadrature points
        function generate_volume_quadrature(obj)
            % Allocate memory space
            % ------------------------------------------------------------------
            nqx = obj.NumberRefTriNodes * obj.NumberVertices;
            obj.CellQuadNodes   = zeros(nqx, 2);
            obj.CellQuadWeights = zeros(nqx, 1);
            % Loop through sub-triangles
            % ------------------------------------------------------------------
            for i=1:obj.NumberVertices
                ind = ((i-1)*obj.NumberRefTriNodes+1):i*obj.NumberRefTriNodes;
                ii = [i,mod(i,obj.NumberVertices)+1];
                vv = [obj.PolyVertices(ii,:);[0,0]]; v0 = vv(1,:);
                J = get_simplex_jacobian(2, vv);
                detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);  svol = detJ/2;
                obj.CellQuadWeights(ind) = obj.RefTriWeights*svol;
                obj.CellQuadNodes(ind,:) = ones(obj.NumberRefTriNodes,1)*v0 + (J*obj.RefTriNodes')';
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Build all surface quadrature points
        function generate_surface_quadrature(obj)
            % Allocate memory space
            % ------------------------------------------------------------------
            obj.FaceQuadNodes   = cell(obj.NumberVertices,1);
            obj.FaceQuadWeights = cell(obj.NumberVertices,1);
            obj.FaceBasisValues = cell(obj.NumberVertices,1);
            obj.FaceBasisGrads  = cell(obj.NumberVertices,1);
            % Loop through faces
            % ------------------------------------------------------------------
            for f=1:obj.NumberVertices
                % Build face quadrature
                fv = obj.PolyFaceNodes{f};
                fverts = obj.PolyVertices(fv,:);
                v0 = fverts(1,:); df = diff(fverts); len = norm(df);
                obj.FaceQuadWeights{f} = len*obj.RefLineWeights;
                obj.FaceQuadNodes{f} = ones(obj.NumberLineNodes,1)*v0 + obj.RefLineNodes*df;
            end     
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate reference polygon basis functions/gradients
        function calculate_refpoly_bfgs(obj)
            % Calculate volume basis functions/gradients
            [bv, gv] = obj.BasisFunc(obj.PolyVertices,obj.CellQuadNodes,obj.PolyFaceNodes,obj.BasisFunctionOrder,obj.NumberVertices);
            obj.CellBasisValues = bv;
            obj.CellBasisGrads  = gv;
            % Calculate surface basis functions/gradients
            for f=1:length(obj.PolyFaceNodes)
                [bs,gs] = obj.BasisFunc(obj.PolyVertices,obj.FaceQuadNodes{f},obj.PolyFaceNodes,obj.BasisFunctionOrder,obj.NumberVertices);
                obj.FaceBasisValues{f} = bs;
                obj.FaceBasisGrads{f}  = gs;
            end
        end
    end
end