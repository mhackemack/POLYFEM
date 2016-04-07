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
            obj.NumberVertices     = varargin{1};
            obj.BasisFunctionName  = varargin{2};
            obj.BasisFunctionOrder = varargin{3};
            obj.QuadratureOrder    = varargin{4};
            % Build Unit Line
            % ------------------------------------------------------------------
            [obj.RefLineNodes, obj.RefLineWeights] = get_legendre_gauss_quad(obj.QuadratureOrder);
            % Build Reference Triangle
            % ------------------------------------------------------------------
            [obj.RefTriNodes, obj.RefTriWeights] = TriGaussPoints(obj.QuadratureOrder);
            obj.NumberRefTriNodes = length(obj.RefTriWeights);
            % Build Reference Polygon
            % ------------------------------------------------------------------
            [obj.PolyVertices,obj.PolyFaceNodes] = RegularPolygon(obj.NumberVertices,1/2);
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
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                          Public Accessor Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        
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
            obj.CellBasisValues = zeros(nqx, obj.NumberVertices);
            obj.CellBasisGrads  = zeros(obj.NumberVertices, 2, nqx);
            % Loop through sub-triangles
            % ------------------------------------------------------------------
            for i=1:obj.NumberVertices
                ii = [i,mod(i,obj.NumberVertices)+1];
                vv = [obj.PolyVertices(ii,:);[0,0]];
                J = get_simplex_jacobian(2, vv);
                detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
                invJ = [J(2,2),-J(1,2);-J(2,1),J(1,1)]/detJ;
                
            end
            % Build basis functions/gradients
            % ------------------------------------------------------------------
            
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
                fv = obj.PolyFaceNodes{f};
                fverts = obj.PolyVertices(fv,:);
                
            end     
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate reference polygon basis functions/gradients
        function calculate_ref_bfgs(obj)
            % Calculate volume basis functions/gradients
            [bv, gv] = obj.BasisFunc();
        end
    end
end