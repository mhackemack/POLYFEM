%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Refine 1D Mesh
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to refine a 1D mesh.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_1D_mesh(mesh)
% Quick Error Checking
% ------------------------------------------------------------------------------
if ~isa(mesh, 'AMRGeometry'), error('Requires AMRGeometry object.'); end
% ------------------------------------------------------------------------------
% Get preliminary information
% ------------------------------------------------------------------------------
num_old_cells = mesh.TotalCells;
num_new_cells = sum(mesh.CellRefinementFlag);
num_new_faces = num_new_cells;
num_new_verts = num_new_cells;
mesh.allocate_more_memory(num_new_verts, num_new_cells, num_new_faces);
% ------------------------------------------------------------------------------
% Loop through cells and refine cell-by-cell
% ------------------------------------------------------------------------------
cc = 0;
for c=1:num_old_cells
    if mesh.CellRefinementFlag(c)
        cc = cc + 1;
        refine_individual_cell(mesh, c, num_old_cells + cc);
    end
end
% ------------------------------------------------------------------------------
% Update Final geometry information
% ------------------------------------------------------------------------------
mesh.update_geometry_info_after_modifications();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_individual_cell( mesh, c, num )
fvnum = num + 1;
m = mesh.MatID(c);
cv1  = mesh.CellVerts{c}(1);  cv2  = mesh.CellVerts{c}(2);
v1   = mesh.Vertices(cv1,:);  v2   = mesh.Vertices(cv2,:);
f1   = mesh.CellFaces{c}(1);  f2   = mesh.CellFaces{c}(2);
fc   = (v1+v2)/2;
% Face Cells
f1c1 = mesh.FaceCells(f1,1); f1c2 = mesh.FaceCells(f1,2);
f2c1 = mesh.FaceCells(f2,1); f2c2 = mesh.FaceCells(f2,2);
% ------------------------------------------------------------------------------
% Make Modifications
% ------------------------------------------------------------------------------
mesh.Vertices(fvnum) = fc;
mesh.MatID(num) = m;
mesh.FaceVerts{fvnum} = fvnum;
mesh.FaceID(fvnum) = 0;
mesh.FaceCells(fvnum,:) = [c,num];
mesh.CellVerts{c} = [cv1,fvnum];   cv1 = mesh.Vertices(mesh.CellVerts{c});
mesh.CellVerts{num} = [fvnum,cv2]; cv2 = mesh.Vertices(mesh.CellVerts{num});
mesh.CellCenter(c,:) = mean(mesh.Vertices(mesh.CellVerts{c}));
mesh.CellCenter(num,:) = mean(mesh.Vertices(mesh.CellVerts{num}));
if cv1(2) < cv1(1), fliplr(mesh.CellVerts{c}); end
if cv2(2) < cv2(1), fliplr(mesh.CellVerts{num}); end
mesh.CellFaces{c} = [f1,fvnum];
mesh.CellFaces{num} = [fvnum,f2];
if f2c1 == c
    mesh.FaceCells(f2,1) = num;
elseif f2c2 == c
    mesh.FaceCells(f2,2) = num;
end
mesh.PreviousCell(c) = c;
mesh.PreviousCell(num) = c;
mesh.CellRefinementLevel(c) = mesh.CellRefinementLevel(c) + 1;
mesh.CellRefinementLevel(num) = mesh.CellRefinementLevel(c);
mesh.CellRefinementTop(num) = mesh.CellRefinementTop(c);
% ------------------------------------------------------------------------------
% Update Refinement Tree - this is actually mostly general for all mesh types
% ------------------------------------------------------------------------------
ctop = mesh.CellRefinementTop(c);
thier = mesh.CellRefinementTreeHierarchy{c}; nthier = length(thier); chier = cell(nthier+1, 1);
ttree = mesh.CellRefinementTree; tncells = {c,num}; tt = ttree; chier{1} = ttree;
% Build new hierarchy tree
for i=1:nthier-1
    ii = i + 1;
    chier{ii,1} = tt{thier(i)};
    tt = tt{thier(i)};
end
chier{end,1} = tncells; tt = tncells;
for i=nthier:-1:1
    ii = thier(i);
    tnew = chier{i};
    tnew{ii} = tt;
    tt = tnew;
end
mesh.CellRefinementTree = tt;
mesh.CellRefinementTreeHierarchy{c} = [thier,1];
mesh.CellRefinementTreeHierarchy{num} = [thier,2];
mesh.CellRefinementTreeHierarchyLevel(ctop) = length(mesh.CellRefinementTreeHierarchy{c}) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%