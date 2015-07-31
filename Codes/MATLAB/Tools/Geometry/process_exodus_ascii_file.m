%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Exodus Mesh File Input
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the Maximum-Entropy basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = process_exodus_ascii_file( varargin )
obj.Dimension = varargin{1};
s_name = varargin{2};
if obj.Dimension == 2
    obj = process_2D_exodus_ascii_file( s_name );
elseif obj.Dimension == 3
    obj = process_3D_exodus_ascii_file( s_name );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Dimension specific routines
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = process_2D_exodus_ascii_file( s_name )
h_name = strcat(s_name, '_header.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = process_3D_exodus_ascii_file( s_name )
h_name = strcat(s_name, '_header.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = allocate_memory_space( obj )
% Vertex Array Structures
obj.Vertices = zeros( obj.TotalVertices, obj.Dimension );
% Cell Array Structures
obj.CellVerts = cell( obj.TotalCells, 1);
obj.CellFaces = cell( obj.TotalCells, 1);
obj.CellVolume = zeros( obj.TotalCells, obj.Dimension );
obj.CellSurfaceArea = zeros( obj.TotalCells, 1 );
% Face Array Structures

% 3D Array Structures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%