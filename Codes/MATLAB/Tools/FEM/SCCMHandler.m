%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Schwarz-Christoffel Conformal Mapper
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:         1) RefPolygonDiskMapper class
%                   2) Global polygon vertices [Nv,2]
%                   3) Global polygon face vertices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef SCCMHandler < handle
    % General SCCMHandler Variables
    properties (Access = public)
        NumberVertices
        Vertices
        FaceVerts
        BasisFunctionName
        BasisFunctionOrder
        PolynomialOrder
    end
    % Unit Disk Variables
    properties (Access = public)
        DiskMap
        DiskVolumeQuadNodes
        DiskVolumeQuadWeights
        DiskVolumeBasisFunctions
        DiskVolumeBasisGradients
        DiskSurfaceQuadNodes
        DiskSurfaceQuadWeights
        DiskSurfaceBasisFunctions
        DiskSurfaceBasisGradients
    end
    % 
    properties (Access = public)
        PolyVolumeQuadNodes
        PolyVolumeQuadWeights
        PolyVolumeBasisFunctions
        PolyVolumeBasisGradients
        PolySurfaceQuadNodes
        PolySurfaceQuadWeights
        PolySurfaceBasisFunctions
        PolySurfaceBasisGradients
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                               Constructor
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        function obj = SCCMHandler (varargin)
            % Retrieve polygon data
            % ------------------------------------------------------------------
            obj.Vertices = varargin{2};
            obj.NumberVertices = size(obj.Vertices,1);
            obj.FaceVerts = varargin{3};
            % Retrieve basis function data
            % ------------------------------------------------------------------
            obj.BasisFunctionName  = varargin{1}.BasisFunctionName;
            obj.BasisFunctionOrder = varargin{1}.BasisFunctionOrder;
            obj.PolynomialOrder    = varargin{1}.PolynomialOrder;
            % Retrieve unit disk data
            % ------------------------------------------------------------------
            obj.DiskMap = varargin{1}.DiskMap;
            
            % 
            % ------------------------------------------------------------------
            
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
        
    end
end