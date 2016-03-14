%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set Linear Transport Solution Face Flags
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = set_quad_sol_face_flags(mesh)
% Get min/max bounds
minX = mesh.minX; maxX = mesh.maxX;
minY = mesh.minY; maxY = mesh.maxY;
% Set boundary flags
mesh.set_face_flag_on_surface(2,[minX,minY;minX,maxY]);
mesh.set_face_flag_on_surface(3,[minX,minY;maxX,minY]);
mesh.set_face_flag_on_surface(3,[minX,maxY;maxX,maxY]);