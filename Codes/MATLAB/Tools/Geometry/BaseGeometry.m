%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Base Geometry Superclass
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
classdef BaseGeometry < handle
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
end