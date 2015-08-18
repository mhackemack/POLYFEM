function write_output_to_vtk_rev2 (filename, data, mesh, DoF, sol, sol_name)
% Quick Error Checking
% --------------------
if nargin < 4, error('Bad input parameters.'); end
if nargin < 5
    sol = [];
    sol_name = [];
end

% Get Basic Information
% ---------------------
dim = mesh.Dimension;
refbool = data.problem.refineMesh;
[y, m, d, h, mi, ~] = datevec(now);

if refbool
    reflvl = mesh.MeshRefinementLevel;
    filename = sprintf('%s_cyc%2.2d',filename,reflvl);
end

% Set VTK filenames
str_mesh = sprintf('%s_mesh.vtk',filename);
str_sol = sprintf('%s_solution.vtk',filename);

% Switch function call based on problem dimensionality
if dim == 1
    write_1D_output_to_vtk(fid1,fid2,mesh,DoF,sol,sol_name,str_mesh,str_sol);
elseif dim == 2
    write_2D_output_to_vtk(fid1,fid2,mesh,DoF,sol,sol_name,str_mesh,str_sol);
elseif dim==3
    write_3D_output_to_vtk(fid1,fid2,mesh,DoF,sol,sol_name,str_mesh,str_sol);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxilliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_1D_output_to_vtk(mesh,DoF,sol,sol_name,str_mesh,str_sol)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_2D_output_to_vtk(mesh,DoF,sol,sol_name,str_mesh,str_sol)

fid2 = generate_solution_file(str_sol);

% Write Mesh Output
% -----------------
% print vertex information
% dim = mesh.Dimension;
% fprintf(fid1,'POINTS %d float\n',mesh.TotalVertices);
% fprintf(fid1,'%f %f %f \n',[mesh.Vertices zeros(mesh.TotalVertices,3-dim)]');
% % print cell information
% ntotdofs = get_total_dof(DoF);
% fprintf(fid1,'CELLS %d %d \n',mesh.TotalCells,mesh.TotalCells+ntotdofs);
% for c=1:mesh.TotalCells
%     cvs = mesh.CellVerts{c};
%     ncvs = length(cvs);
%     fprintf(fid1,' %d ',ncvs);
%     for k=1:ncvs
%         fprintf(fid1,'%d ',cvs(k)-1);
%     end
%     fprintf(fid1,' \n');
% end
% fprintf(fid1,' \n');
% % print cell types
% fprintf(fid1,'CELL_TYPES %d\n',mesh.TotalCells);
% fprintf(fid1,'%d\n',7*ones(mesh.TotalCells,1));

% Write Solution Output
% ---------------------
if isempty(sol), return; end
if ~iscell(sol_name)
    sol_name = {sol_name};
end
if length(sol_name) > 1 && length(sol) > 1 && iscell(sol)
    if length(sol_name) ~= length(sol)
        error('Length of sol structure does not equal length of sol_name structure.')
    end
elseif length(sol_name) == 1
    if ~iscell(sol)
        tsol = sol; clear sol;
        sol{1} = tsol;
    end
    if ~iscell(sol_name)
        tsol_name = sol_name; clear sol_name;
        sol_name{1} = tsol_name;
    end
else
    error('Solution structures are weird...')
end
% print vertex information
fprintf(fid2,'POINTS %d float \n',DoF.TotalDoFs);
fprintf(fid2,'%f %f %f \n',[DoF.NodeLocations zeros(DoF.TotalDoFs,3-dim)]');
% print cell information
fprintf(fid2,'CELLS %d %d \n',DoF.TotalCells,DoF.TotalCells+ntotdofs);
for c=1:DoF.TotalCells
    cvs = DoF.ConnectivityArray{c};
    ncvs = length(cvs);
    fprintf(fid2,' %d ',ncvs);
    for k=1:ncvs
        fprintf(fid2,'%d ',cvs(k)-1);
    end
    fprintf(fid2,' \n');
end
fprintf(fid2,' \n');
% print cell types
fprintf(fid2,'CELL_TYPES %d\n',DoF.TotalCells);
fprintf(fid2,'%d\n',7*ones(DoF.TotalCells,1));
% add solution data
fprintf(fid2,'POINT_DATA %d %d \n',DoF.TotalDoFs);
for s=1:length(sol_name)
    s_name = sol_name{s};
    ssol = sol{s};
    fprintf(fid2,'SCALARS %s double \n',s_name);
    fprintf(fid2,'LOOKUP_TABLE   default \n');
    fprintf(fid2,'%f\n',ssol);
    fprintf(fid2, ' \n');
end

fclose(fid2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_3D_output_to_vtk(mesh,DoF,sol,sol_name,str_mesh,str_sol)

fid1 = generate_mesh_file(str_mesh);
fid2 = generate_solution_file(str_sol);

% Write Mesh Output
% -----------------
% Get some geometry info
v_dofs = []; t_dofs = [];
num_mesh_cells = 0; num_cells = 0;
for c=1:mesh.TotalCells
    cdofs = DoF.ConnectivityArray{c};
    cfaces = mesh.CellFaces{c};
    cmean = mean(DoF.NodeLocations(cdofs,:));
    for ff=1:length(cfaces)
        fdofs = DoF.CellFaceNodes{c}{ff};
        fmean = mean(DoF.NodeLocations(fdofs,:));
        num_mesh_cells = num_mesh_cells + 1;
        v_dofs = [v_dofs;DoF.NodeLocations(fdofs,:)];
        for i=1:length(fdofs)
            if i==length(fdofs)
                ii = [i,1];
            else
                ii = [i,i+1];
            end
            num_cells = num_cells + 1;
            t_dofs = [t_dofs;DoF.NodeLocations(fdofs(ii),:);fmean;cmean];
        end
    end
end
num_dofs = size(v_dofs,1);
% Print Points Information
fprintf(fid1,'POINTS %d float \n',num_dofs);
fprintf(fid1,'%f %f %f \n',v_dofs');
% Print Cell Information
fprintf(fid1,'CELLS %d %d \n',num_mesh_cells,num_mesh_cells+num_dofs);
n=0;
for c=1:mesh.TotalCells
    cfaces = mesh.CellFaces{c};
    for ff=1:length(cfaces)
        fdofs = DoF.CellFaceNodes{c}{ff};
        fprintf(fid1,' %d ',length(fdofs));
        for i=1:length(fdofs)
            fprintf(fid1,' %d ',n);
            n = n + 1;
        end
        fprintf(fid1,' \n');
    end
end
fprintf(fid1,' \n');
% print cell types
fprintf(fid1,'CELL_TYPES %d\n',num_mesh_cells);
fprintf(fid1,'%d\n',7*ones(num_mesh_cells,1));

% Write Solution Output
% ---------------------
if isempty(sol), return; end
if ~iscell(sol_name)
    sol_name = {sol_name};
end
if length(sol_name) > 1 && length(sol) > 1 && iscell(sol)
    if length(sol_name) ~= length(sol)
        error('Length of sol structure does not equal length of sol_name structure.')
    end
elseif length(sol_name) == 1
    if ~iscell(sol)
        sol = {sol};
    end
    if ~iscell(sol_name)
        tsol_name = sol_name; clear sol_name;
        sol_name{1} = tsol_name;
    end
else
    error('Solution structures are weird...')
end

% Get some geometry info
% t_dofs = [];
% num_cells = 0; 
% for c=1:mesh.TotalCells
%     cdofs = DoF.ConnectivityArray{c};
%     cmean = mean(DoF.NodeLocations(cdofs,:));
%     cfaces = mesh.CellFaces{c};
%     for ff=1:length(cfaces)
%         fdofs = DoF.CellFaceNodes{c}{ff};
%         fmean = mean(DoF.NodeLocations(fdofs,:));
%         for i=1:length(fdofs)
%             if i==length(fdofs)
%                 ii = [i,1];
%             else
%                 ii = [i,i+1];
%             end
%             num_cells = num_cells + 1;
%             t_dofs = [t_dofs;DoF.NodeLocations(fdofs(ii),:);fmean;cmean];
%         end
%     end
% end
num_dofs = size(t_dofs,1);
% Print Points Information
fprintf(fid2,'POINTS %d float \n',num_dofs);
fprintf(fid2,'%f %f %f \n',t_dofs');
% Print Cell Information
fprintf(fid2,'CELLS %d %d \n',num_cells,num_cells*5);
n = 0;
for c=1:num_cells
    fprintf(fid2,' %d ',4);
    for i=1:4
%         n = n + 1;
        fprintf(fid2,' %d ',n);
        n = n + 1;
    end
    fprintf(fid2,' \n');
end
fprintf(fid2,' \n');
% print cell types
fprintf(fid2,'CELL_TYPES %d\n',num_cells);
fprintf(fid2,'%d\n',10*ones(num_cells,1));
% add solution data
fprintf(fid2,'POINT_DATA %d %d \n',num_dofs);
n = 0;
for s=1:length(sol_name)
    s_name = sol_name{s};
    tsol = zeros(num_dofs, 1);
    fprintf(fid2,'SCALARS %s double \n',s_name);
    fprintf(fid2,'LOOKUP_TABLE   default \n');
    for c=1:mesh.TotalCells
        cdofs = DoF.ConnectivityArray{c};
        cmean = mean(sol{s}(cdofs));
        cfaces = mesh.CellFaces{c};
        for ff=1:length(cfaces)
            fdofs = DoF.CellFaceNodes{c}{ff};
            fmean = mean(sol{s}(fdofs));
            for i=1:length(fdofs)
                if i==length(fdofs)
                    ii = [i,1];
                else
                    ii = [i,i+1];
                end 
                td = fdofs(ii);
                tsol(n+1:n+4) = [sol{s}(td);fmean;cmean];
                n = n + 4;
            end
        end
    end
    fprintf(fid2,'%f\n',tsol);
end

fclose(fid1);
fclose(fid2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_total_dof(DoF)
if DoF.FEMType == 1
    out = 0;
    for c=1:DoF.TotalCells
        out = out + length(DoF.ConnectivityArray{c});
    end
else
    out = DoF.TotalDoFs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fid = generate_mesh_file(filename)
[y, m, d, h, mi, ~] = datevec(now);
str_mesh = sprintf('%s_mesh.vtk',filename);
fid = fopen(str_mesh,'w');
% Print Basic Information
fprintf(fid,'# vtk DataFile Version 3.0 \n');
fprintf(fid,'Date: %d/%d/%d   Time: %d:%d\n', m, d, y, h, mi);
fprintf(fid,'ASCII \n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
fprintf(fid,' \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fid = generate_solution_file(filename)
[y, m, d, h, mi, ~] = datevec(now);
str_sol = sprintf('%s_solution.vtk',filename);
fid = fopen(str_sol,'w');
% Print Basic Information
fprintf(fid,'# vtk DataFile Version 3.0 \n');
fprintf(fid,'Date: %d/%d/%d   Time: %d:%d\n', m, d, y, h, mi);
fprintf(fid,'ASCII \n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
fprintf(fid,' \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%