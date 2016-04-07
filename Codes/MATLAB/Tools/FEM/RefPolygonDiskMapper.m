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
%                   4) Polynomial Completeness Order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef RefPolygonDiskMapper < handle
    % General Mapper Variables
    properties (Access = public)
        NumberVertices
        BasisFunctionName
        BasisFunctionOrder
        PolynomialOrder
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
            obj.PolynomialOrder    = varargin{4};
            % Build Unit Line
            % ------------------------------------------------------------------
            [obj.RefLineNodes, obj.RefLineWeights] = get_legendre_gauss_quad(obj.PolynomialOrder);
            % Build Reference Triangle
            % ------------------------------------------------------------------
            [obj.RefTriNodes, obj.RefTriWeights] = TriGaussPoints(obj.PolynomialOrder);
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
        % Calculate reference polygon basis functions/gradients
        function calculate_ref_bfgs(obj)
            % Calculate volume basis functions/gradients
            [bv, gv] = obj.BasisFunc();
        end
    end
end