%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Refine Triangle Mesh
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to refine a 2D triangular mesh.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_triangle_mesh(mesh)
% Quick Error Checking
% ------------------------------------------------------------------------------
if ~isa(mesh, 'AMRGeometry'), error('Requires AMRGeometry object.'); end
% ------------------------------------------------------------------------------
% Sort cells from lowest refinement level - this yields more consistency in
% refinement routines (I think...).
% ------------------------------------------------------------------------------
% No limit on AMR Irregularity
if isinf(mesh.MaxIrregularity)
    new_cells = (1:mesh.TotalCells)';
    new_cells = new_cells(mesh.CellRefinementFlag);
    current_lvls = mesh.CellRefinementLevel(new_cells);
    [~, ind] = sort(current_lvls);
    new_cells = new_cells(ind); num_new_cells = length(new_cells);
% Impose AMR Irregularity limits
else
    ref_cells = (1:mesh.TotalCells)'; ref_cells(~mesh.CellRefinementFlag) = [];
    next_lvls = mesh.CellRefinementLevel;
    next_lvls(ref_cells) = next_lvls(ref_cells) + 1;
    [next_lvls, new_cells] = mesh.check_cell_refinement_differences(next_lvls, ref_cells);
    [~, ind] = sort(next_lvls(new_cells));
    new_cells = new_cells(ind); num_new_cells = length(new_cells);
end
disp(['   -> Number of Refinement Flags: ',num2str(num_new_cells)])
% ------------------------------------------------------------------------------
% Loop through cells and refine cell-by-cell
% ------------------------------------------------------------------------------
c_count = mesh.TotalCells; f_count = mesh.TotalFaces; v_count = mesh.TotalVertices;
mesh.allocate_more_memory(0,num_new_cells*3,0);
mesh.CellRefinedLastCycle = false(mesh.TotalCells,1);
rev_str = [];
for c=1:length(new_cells)
    msg = sprintf('      -> Refining Cell: %d of %d',c,num_new_cells);
    fprintf([rev_str,msg]);
    rev_str = repmat(sprintf('\b'), 1, length(msg));
    tcell = new_cells(c);
    % Determine new cell/face/vert counts
    ncells = 3; nfaces = 3; nverts = 0;
    for ff=1:length(mesh.RefCellHigherLvls{tcell})
        if ~mesh.RefCellHigherLvls{tcell}(ff)
            nfaces = nfaces + 1;
            nverts = nverts + 1;
        end
    end
    % Allocate more array memory for mesh
    mesh.allocate_more_memory(nverts,0,nfaces);
    % Refine individual cell
    cnums = c_count + (1:ncells);
    fnums = f_count + (1:nfaces);
    vnums = v_count + (1:nverts);
    refine_individual_cell(mesh, tcell, cnums, fnums, vnums);
    % Update Counts
    c_count = c_count + ncells;
    f_count = f_count + nfaces;
    v_count = v_count + nverts;
    mesh.TotalCells = mesh.TotalCells + ncells;
    mesh.TotalFaces = mesh.TotalFaces + nfaces;
    mesh.TotalVertices = mesh.TotalVertices + nverts;
end
% ------------------------------------------------------------------------------
% Update Final geometry information
% ------------------------------------------------------------------------------
fprintf(rev_str);
mesh.update_geometry_info_after_modifications();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine Individual Cell
%   Variable Listing:
%   1) mesh  - Reference to AMRGeometry object
%   2) c     - Cell to be refined
%   3) cnums - Cell IDs for new cells to be created
%   4) fnums - Face IDs for new faces to be created
%   5) vnums - Vertex IDs for new vertices to be created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refine_individual_cell( mesh, c, cnums, fnums, vnums )
% Loop through macro faces and retrieve information
% ------------------------------------------------------------------------------
ncells = length(cnums); nfaces = length(fnums); nverts = length(vnums);
current_lvl = mesh.CellRefinementLevel(c);
next_lvl = current_lvl + 1;
cell_macro_faces = mesh.RefCellFaces{c};
cell_macro_face_verts = mesh.RefCellFaceVerts{c};
cell_macro_face_nums = mesh.RefCellFaceNumbering{c};
cell_macro_cells = mesh.RefCellFaceCells{c};
cell_corner_verts = mesh.RefCellCornerVerts{c};
cell_mid_face_verts = mesh.RefCellMidFaceVerts{c};
cell_center = mean(mesh.Vertices(cell_corner_verts,:));
has_higher_lvl_cells = mesh.RefCellHigherLvls{c};
% Some vertex arrays
new_verts = zeros(nverts, mesh.Dimension);
needed_vert_inds = zeros(6,1);
needed_vert_inds(1:3) = cell_corner_verts;
needed_verts = zeros(6,mesh.Dimension);
needed_verts(1:3,:) = mesh.Vertices(cell_corner_verts,:);
% Some face arrays
cell_boundary_faces = false(length(mesh.RefCellFaces{c}),1);
new_face_vert_inds = zeros(nfaces, 2);
new_face_ids = zeros(nfaces, 1);
new_old_face_vert_inds = [];
new_old_face_combo = zeros(3,2);
% Loop through macro faces
nv_count = 1;
for ff=1:length(cell_macro_faces)
    of = mesh.RefCellFaces{c}{ff};
    fvind = [ff,mod(ff,length(cell_macro_faces))+1];
    fcorn_verts = cell_corner_verts(fvind);
    tcs = mesh.RefCellFaceCells{c}{ff};
    new_face_ids(ff) = 0;
    new_face_vert_inds(ff,:) = ff;
    % Boundary Face
    if isempty(tcs)
        cell_boundary_faces(ff) = true;
        new_verts(nv_count,:) = mean(mesh.Vertices(fcorn_verts,:));
        needed_verts(3+ff,:) = mean(mesh.Vertices(fcorn_verts,:));
        needed_vert_inds(3+ff) = vnums(nv_count);
        new_face_vert_inds(3+nv_count,1) = 3+ff;
        new_face_vert_inds(3+nv_count,2) = ff;
        new_face_ids(3+nv_count) = mesh.FaceID(of);
        new_old_face_vert_inds = [new_old_face_vert_inds;ff,fcorn_verts(1),vnums(nv_count)];
        new_old_face_combo(ff,:) = [cell_macro_faces{ff},fnums(3+nv_count)];
        nv_count = nv_count + 1;
    % Face with refinement at same or lower level on other side
    elseif ~has_higher_lvl_cells(ff)
        new_verts(nv_count,:) = mean(mesh.Vertices(fcorn_verts,:));
        needed_verts(3+ff,:) = mean(mesh.Vertices(fcorn_verts,:));
        needed_vert_inds(3+ff) = vnums(nv_count);
        new_face_vert_inds(3+nv_count,1) = 3+ff;
        new_face_vert_inds(3+nv_count,2) = ff;
        new_face_ids(3+nv_count) = mesh.FaceID(of);
        new_old_face_vert_inds = [new_old_face_vert_inds;ff,fcorn_verts(1),vnums(nv_count)];
        new_old_face_combo(ff,:) = [cell_macro_faces{ff},fnums(3+nv_count)];
        nv_count = nv_count + 1;
    % Face with higher refinement on other side
    elseif has_higher_lvl_cells(ff)
        needed_verts(3+ff,:) = mesh.Vertices(cell_mid_face_verts(ff),:);
        needed_vert_inds(3+ff) = cell_mid_face_verts(ff);
    end
end
needed_vert_inds = needed_vert_inds';
new_cell_corner_vert_inds = get_new_cell_corner_vertex_indices();
new_cell_macro_faces = get_new_cell_ext_faces();
new_cell_int_faces = get_new_cell_int_faces();
new_face_inds = get_new_face_indices();
rep_cell_nums = get_replacement_cell_nums();
rep_face_nums = get_replacement_face_nums();
% ------------------------------------------------------------------------------
% Make Modifications for Particular Cell/Faces
% ------------------------------------------------------------------------------
tcnums = [c,cnums];
mesh.MatID(tcnums) = mesh.MatID(c);
% Clear Ref Cell Arrays
for cc=1:length(tcnums)
    mesh.RefCellFaces{tcnums(cc)} = cell(3,1);
    mesh.RefCellFaceCells{tcnums(cc)} = cell(3,1);
    mesh.RefCellFaceVerts{tcnums(cc)} = cell(3,1);
    mesh.RefCellMidFaceVerts{tcnums(cc)} = zeros(1,3);
    mesh.RefCellHigherLvls{tcnums(cc)} = false(1,3);
    mesh.RefCellFaceNumbering{tcnums(cc)} = cell(3,1);
end
mesh.Vertices(vnums,:) = new_verts;
mesh.FaceID(fnums) = new_face_ids;
% Modify newly created faces
for ff=1:nfaces
    f = fnums(ff);
    if ff <= 3
        tfv = needed_vert_inds(new_face_inds(ff,:));
        mesh.FaceCells(f,:) = [tcnums(ff),tcnums(end)];
    else
        iif = new_face_vert_inds(ff,2);
        iiff = mod(iif,3)+1;
        tcmf = new_face_vert_inds(ff,1);
        tfv = needed_vert_inds(new_face_inds(tcmf,:));
        mesh.FaceCells(f,:) = zeros(1,2);
        mesh.FaceCells(f,1) = tcnums(iiff);
        if ~cell_boundary_faces(iif)
            mesh.FaceCells(f,2) = cell_macro_cells{iif};
        end
        oface = cell_macro_faces{iif};
        mesh.FaceCells(oface,1) = tcnums(iif);
        if mesh.FaceID(oface) == 0
            mesh.FaceCells(oface,2) = cell_macro_cells{iif};
        end
        ntfv = needed_vert_inds([iif,3+iif]);
        mesh.FaceVerts{oface} = ntfv;
        mesh.FaceCenter(oface,:) = mean(mesh.Vertices(ntfv,:));
    end
    mesh.FaceVerts{f} = tfv;
    mesh.FaceCenter(f,:) = mean(mesh.Vertices(tfv,:));
end
% Loop through cells
for cc=1:length(tcnums)
    ncif = new_cell_int_faces{cc};
    mesh.RefCellCornerVerts{tcnums(cc)} = needed_vert_inds(new_cell_corner_vert_inds(cc,:));
    % Modify Interior Faces
    % ---------------------
    for i=1:size(ncif,1)
        tncif = ncif(i,:);
        mesh.RefCellFaces{tcnums(cc)}{tncif(1)} = fnums(tncif(2));
        mesh.RefCellFaceCells{tcnums(cc)}{tncif(1)} = tcnums(tncif(3));
        mesh.RefCellMidFaceVerts{tcnums(cc)}(tncif(1)) = 0;
        mesh.RefCellFaceNumbering{tcnums(cc)}{tncif(1)} = tncif(4);
        % Inside Cell
        if cc==4
            mesh.RefCellFaceVerts{tcnums(cc)}{tncif(1)} = fliplr(mesh.FaceVerts{fnums(tncif(2))});
        % Outer Cells
        else
            mesh.RefCellFaceVerts{tcnums(cc)}{tncif(1)} = mesh.FaceVerts{fnums(tncif(2))};
        end
    end
    % Modify Exterior Faces
    % ---------------------
    if cc ~= 4
        f1 = new_cell_macro_faces(cc,1); f2 = new_cell_macro_faces(cc,2);
        % First Exterior Face
        if ~has_higher_lvl_cells(f1)
            mesh.RefCellFaces{tcnums(cc)}{1} = new_old_face_combo(f1,1);
            mesh.RefCellFaceCells{tcnums(cc)}{1} = cell_macro_cells{f1};
            mesh.RefCellFaceVerts{tcnums(cc)}{1} = needed_vert_inds([f1,3+f1]);
            mesh.RefCellFaceNumbering{tcnums(cc)}{1} = cell_macro_face_nums{f1};
        else
            tcv = cell_macro_face_verts{f1};
            tcf = cell_macro_faces{f1};
            tcc = cell_macro_cells{f1};
            mcfv = cell_mid_face_verts(f1);
            tcfn = cell_macro_face_nums{f1};
            for i=1:length(tcv)
                if mcfv == tcv(i)
                    if mcfv == tcv(i)
                        mesh.RefCellFaceVerts{tcnums(cc)}{1} = tcv(1:i);
                        mesh.RefCellFaces{tcnums(cc)}{1} = tcf(1:i-1);
                        mesh.RefCellFaceCells{tcnums(cc)}{1} = tcc(1:i-1);
                        mesh.RefCellFaceNumbering{tcnums(cc)}{1} = tcfn(1:i-1);
                        if length(mesh.RefCellFaceCells{tcnums(cc)}{1}) == 1
                            mesh.RefCellMidFaceVerts{tcnums(cc)}(1) = 0;
                            mesh.RefCellHigherLvls{tcnums(cc)}(1) = false;
                        else
                            mesh.RefCellHigherLvls{tcnums(cc)}(1) = true;
                            rcfv = mesh.RefCellFaceVerts{tcnums(cc)}{1};
                            rccv = mesh.RefCellCornerVerts{tcnums(cc)};
                            mrccv = mean(mesh.Vertices(rccv(1:2),:));
                            for ii=1:length(rcfv)
                                if norm(mrccv-mesh.Vertices(rcfv(ii),:)) < 1e-13
                                    mesh.RefCellMidFaceVerts{tcnums(cc)}(1) = rcfv(ii);
                                    break
                                end
                            end
                        end
                        break
                    end
                end
            end
        end
        % Second Exterior Face
        if ~has_higher_lvl_cells(f2)
            mesh.RefCellFaces{tcnums(cc)}{3} = new_old_face_combo(f2,2);
            mesh.RefCellFaceCells{tcnums(cc)}{3} = cell_macro_cells{f2};
            mesh.RefCellFaceVerts{tcnums(cc)}{3} = [needed_vert_inds(3+f2),needed_vert_inds(cc)];
            mesh.RefCellFaceNumbering{tcnums(cc)}{3} = cell_macro_face_nums{f2};
        else
            tcv = cell_macro_face_verts{f2};
            tcf = cell_macro_faces{f2};
            tcc = cell_macro_cells{f2};
            mcfv = cell_mid_face_verts(f2);
            tcfn = cell_macro_face_nums{f2};
            for i=1:length(tcv)
                if mcfv == tcv(i)
                    mesh.RefCellFaceVerts{tcnums(cc)}{3} = tcv(i:end);
                    mesh.RefCellFaces{tcnums(cc)}{3} = tcf(i:end);
                    mesh.RefCellFaceCells{tcnums(cc)}{3} = tcc(i:end);
                    mesh.RefCellFaceNumbering{tcnums(cc)}{3} = tcfn(i:end);
                    if length(mesh.RefCellFaceCells{tcnums(cc)}{3}) == 1
                        mesh.RefCellMidFaceVerts{tcnums(cc)}(3) = 0;
                        mesh.RefCellHigherLvls{tcnums(cc)}(3) = false;
                    else
                        mesh.RefCellHigherLvls{tcnums(cc)}(3) = true;
                        rcfv = mesh.RefCellFaceVerts{tcnums(cc)}{3};
                        rccv = mesh.RefCellCornerVerts{tcnums(cc)};
                        mrccv = mean(mesh.Vertices(rccv([1,3]),:));
                        for ii=1:length(rcfv)
                            if norm(mrccv-mesh.Vertices(rcfv(ii),:)) < 1e-13
                                mesh.RefCellMidFaceVerts{tcnums(cc)}(3) = rcfv(ii);
                                break
                            end
                        end
                    end
                    break
                end
            end
        end
    end
    % Accumulate Cell Faces/Vertices
    cfaces = []; cverts = [];
    for i=1:3
        cfaces = [cfaces,mesh.RefCellFaces{tcnums(cc)}{i}];
        cverts = [cverts,mesh.RefCellFaceVerts{tcnums(cc)}{i}];
    end
    mesh.CellFaces{tcnums(cc)} = unique(cfaces, 'stable');
    mesh.CellVerts{tcnums(cc)} = unique(cverts, 'stable');
    mesh.CellCenter(tcnums(cc),:) = mean(mesh.Vertices(mesh.CellVerts{tcnums(cc)},:));
end

% ------------------------------------------------------------------------------
% Make Modifications for Cell/Faces Neighbors
% ------------------------------------------------------------------------------
for ff=1:3
    if cell_boundary_faces(ff), continue; end
    mc_faces = cell_macro_faces{ff};
    mc_cells = cell_macro_cells{ff};
    mcf_nums = cell_macro_face_nums{ff};
    mc_verts = cell_macro_face_verts{ff};
    % Face with refinement at same or lower level on other side
    if ~has_higher_lvl_cells(ff)
        ocverts = mesh.RefCellFaceVerts{mc_cells}{mcf_nums};
        ocfaces = mesh.RefCellFaces{mc_cells}{mcf_nums};
        occells = mesh.RefCellFaceCells{mc_cells}{mcf_nums};
        oclvl   = mesh.CellRefinementLevel(mc_cells);
        ocnum   = mesh.RefCellFaceNumbering{mc_cells}{mcf_nums};
        mesh.RefCellHigherLvls{mc_cells}(mcf_nums) = true;
        if oclvl == current_lvl
            mesh.RefCellMidFaceVerts{mc_cells}(mcf_nums) = needed_vert_inds(3+ff);
        end
        tc = []; tf = []; tv = ocverts(1); tnum = [];
        for i=1:length(occells)
            if occells(i) == c
                tv = [tv,needed_vert_inds(3+ff),ocverts(i+1:end)];
                tf = [tf,new_old_face_combo(ff,[2,1]),ocfaces(i+1:end)];
                temp_cells = tcnums(rep_cell_nums(ff,:));
                tc = [tc,temp_cells,occells(i+1:end)];
                tnum = [tnum,rep_face_nums(ff,:),ocnum(i+1:end)];
                break
            else
                tf = [tf,ocfaces(i)];
                tv = [tv,ocverts(i+1)];
                tc = [tc,occells(i)];
                tnum = [tnum,ocnum(i)];
            end
        end
        mesh.RefCellFaceVerts{mc_cells}{mcf_nums} = tv;
        mesh.RefCellFaces{mc_cells}{mcf_nums} = tf;
        mesh.RefCellFaceCells{mc_cells}{mcf_nums} = tc;
        mesh.RefCellFaceNumbering{mc_cells}{mcf_nums} = tnum;
        cfaces = []; cverts = [];
        for i=1:3
            cfaces = [cfaces,mesh.RefCellFaces{mc_cells}{i}];
            cverts = [cverts,mesh.RefCellFaceVerts{mc_cells}{i}];
        end
        mesh.CellFaces{mc_cells} = unique(cfaces, 'stable');
        mesh.CellVerts{mc_cells} = unique(cverts, 'stable');
    % Face with higher refinement on other side
    else
        mcfv = cell_mid_face_verts(ff);
        for i=1:length(mc_verts)
            if mc_verts(i) == mcfv
                mcfv_ind = i;
                break;
            end
        end
        % Loop through cells along face
        for i=1:length(mc_cells)
            tcell = mc_cells(i);
            tface = mc_faces(i);
            tnum = mcf_nums(i);
            if i+1 <= mcfv_ind
                t_ind = 2;
                trcfn = rep_cell_nums(ff,2);
            else
                t_ind = 1;
                trcfn = rep_cell_nums(ff,1);
            end
            mesh.RefCellFaceCells{tcell}{tnum} = tcnums(trcfn);
            mesh.RefCellFaceNumbering{tcell}{tnum} = rep_face_nums(ff,t_ind);
            if mesh.FaceCells(tface,1) == c
                mesh.FaceCells(tface,1) = tcnums(trcfn);
            elseif mesh.FaceCells(tface,2) == c
                mesh.FaceCells(tface,2) = tcnums(trcfn);
            end
        end
    end
end
% ------------------------------------------------------------------------------
% Update Refinement Tree - this is actually mostly general for all mesh types
% ------------------------------------------------------------------------------
mesh.PreviousCell(tcnums) = c;
mesh.CellRefinedLastCycle = false(mesh.TotalCells,1);
mesh.CellRefinementLevel(tcnums) = mesh.CellRefinementLevel(c) + 1;
mesh.CellRefinementTop(tcnums) = mesh.CellRefinementTop(c);
ctop = mesh.CellRefinementTop(c);
thier = mesh.CellRefinementTreeHierarchy{c}; nthier = length(thier); chier = cell(nthier+1, 1);
ttree = mesh.CellRefinementTree; tncells{1} = c; tt = ttree; chier{1} = ttree;
for i=1:ncells, tncells{i+1} = cnums(i); end
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
for i=1:ncells, mesh.CellRefinementTreeHierarchy{cnums(i)} = [thier,i+1]; end
mesh.CellRefinementTreeHierarchyLevel(ctop) = length(mesh.CellRefinementTreeHierarchy{c}) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_new_cell_corner_vertex_indices()
out = [1,4,6;2,5,4;3,6,5;4,5,6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_new_face_indices()
out = [4,6;5,4;6,5;4,2;5,3;6,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_new_cell_int_faces()
out{1} = [2,1,4,3];
out{2} = [2,2,4,1];
out{3} = [2,3,4,2];
out{4} = [1,2,2,2;2,3,3,2;3,1,1,2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_new_cell_ext_faces()
out = [1,3;2,1;3,2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_replacement_cell_nums()
out = [2,1;3,2;1,3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_replacement_face_nums()
out = [3,1;3,1;3,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
