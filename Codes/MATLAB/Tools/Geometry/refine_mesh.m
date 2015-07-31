%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Refine Mesh 
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB script to switch between different refinement
%                   methods.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = refine_mesh(mesh)
dim = mesh.Dimension;
mtype = mesh.OriginalMeshType;
if dim == 1
    mesh = refine_1D_mesh(mesh);
elseif dim == 2
    if strcmp(mtype, 'Quadrilateral')
        mesh = refine_quad_mesh(mesh);
    elseif stcmp(mtype, 'Triangle')
        mesh = refine_triangle_mesh(mesh);
    else
        error('Only 2D quad and triangle CURRENTLY supported for refinement.');
    end
elseif dim == 3
    error('3D mesh refinement not yet written.');
else
    error('How did the program get to here...?')
end