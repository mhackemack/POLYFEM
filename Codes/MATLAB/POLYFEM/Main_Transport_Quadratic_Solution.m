%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Spectral Radius Transport Run Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Project Space
% ------------------------------------------------------------------------------
clc; close all; format long e; clear; fclose all;
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Office');
inp = 'Transport_Quad_Sol';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, ~] = load_user_input();
% Begin user input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BF_names = {'PWLD'};
% BF_names = {'PWLD','WACHSPRESS','MV','MAXENT'};
data.problem.Dimension = 2;
print_err_bool = true;
data.problem.saveVTKSolution = 0;
% End user input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartesian Meshes
% ------------------------------------------------------------------------------
% ncells = 10; L = 1;
% x = linspace(0,L,ncells+1);
% geometry = CartesianGeometry(2,x,x);
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     data.problem.Name = sprintf('SOLUTION_cart_%s_k2',now_name);
%     [data, geometry] = process_input_data(data, geometry);
%     data = cleanup_neutronics_input_data(data, geometry);
%     [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
%     % Print difference in solution
%     if print_err_bool
%         data.problem.Name = sprintf('cart_%s_k2',now_name);
%         outname = ['outputs/',data.problem.Path,'/','ERROR_',data.problem.Name];
%         ye = ExactSol_QuadSolution(DoF.NodeLocations);
%         df = ye - sol.flux{:}; df(abs(df)>1e-12) = 0;
%         write_output_to_vtk_rev2(outname,data,geometry,DoF,df,'flux');
%     end
% end
% clear geometry;
% % Triangle Meshes
% % ------------------------------------------------------------------------------
% ncells = 10; L = 1;
% x = linspace(0,L,ncells+1);
% [x,y]=meshgrid(x,x);
% x=x(:); y=y(:);
% tri = delaunayTriangulation(x,y);
% geometry = GeneralGeometry(2, 'Delaunay', tri);
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     data.problem.Name = sprintf('SOLUTION_tri_%s_k2',now_name);
%     [data, geometry] = process_input_data(data, geometry);
%     data = cleanup_neutronics_input_data(data, geometry);
%     [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
%     % Print difference in solution
%     if print_err_bool
%         data.problem.Name = sprintf('tri_%s_k2',now_name);
%         outname = ['outputs/',data.problem.Path,'/','ERROR_',data.problem.Name];
%         ye = ExactSol_QuadSolution(DoF.NodeLocations);
%         df = ye - sol.flux{:}; df(abs(df)>1e-12) = 0;
%         write_output_to_vtk_rev2(outname,data,geometry,DoF,df,'flux');
%     end
% end
% clear geometry;
% % Shestakov Quad Meshes
% % ------------------------------------------------------------------------------
% gname = 'shestakov_quad_L1_nc7_emb3_a0.2';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     data.problem.Name = sprintf('SOLUTION_shes_quad_%s_k2',now_name);
%     [data, geometry] = process_input_data(data, geometry);
%     data = cleanup_neutronics_input_data(data, geometry);
%     [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
%     % Print difference in solution
%     if print_err_bool
%         data.problem.Name = sprintf('shes_quad_%s_k2',now_name);
%         outname = ['outputs/',data.problem.Path,'/','ERROR_',data.problem.Name];
%         ye = ExactSol_QuadSolution(DoF.NodeLocations);
%         df = ye - sol.flux{:}; df(abs(df)>1e-12) = 0;
%         write_output_to_vtk_rev2(outname,data,geometry,DoF,df,'flux');
%     end
% end
% clear geometry;
% % Shestakov Poly Meshes
% % ------------------------------------------------------------------------------
% gname = 'shestakov_poly_mesh_L1_nc3_a0.15';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     data.problem.Name = sprintf('SOLUTION_shes_poly_%s_k2',now_name);
%     [data, geometry] = process_input_data(data, geometry);
%     data = cleanup_neutronics_input_data(data, geometry);
%     [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
%     % Print difference in solution
%     if print_err_bool
%         data.problem.Name = sprintf('shes_poly_%s_k2',now_name);
%         outname = ['outputs/',data.problem.Path,'/','ERROR_',data.problem.Name];
%         ye = ExactSol_QuadSolution(DoF.NodeLocations);
%         df = ye - sol.flux{:}; df(abs(df)>1e-12) = 0;
%         write_output_to_vtk_rev2(outname,data,geometry,DoF,df,'flux');
%     end
% end
% clear geometry;
% % Smooth Poly Meshes
% % ------------------------------------------------------------------------------
gname = 'smooth_poly_mesh_L1_n8_a0.15';
load(strcat(glob.geom_path,gname,'.mat'));
for b=1:length(BF_names)
    now_name = BF_names{b};
    data.Neutronics.SpatialMethod = BF_names{b};
    data.problem.Name = sprintf('SOLUTION_smooth_poly_%s_k2',now_name);
    [data, geometry] = process_input_data(data, geometry);
    data = cleanup_neutronics_input_data(data, geometry);
    [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
    % Print difference in solution
    if print_err_bool
        data.problem.Name = sprintf('smooth_poly_%s_k2',now_name);
        outname = ['outputs/',data.problem.Path,'/','ERROR_',data.problem.Name];
        ye = ExactSol_QuadSolution(DoF.NodeLocations);
        df = ye - sol.flux{:}; df(abs(df)>1e-12) = 0;
        write_output_to_vtk_rev2(outname,data,geometry,DoF,df,'flux');
    end
end
% clear geometry;
% % Z-Quad Meshes
% % ------------------------------------------------------------------------------
% gname = 'z_mesh_quad_L1_n20_a0.2';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     data.problem.Name = sprintf('SOLUTION_z_quad_%s_k2',now_name);
%     [data, geometry] = process_input_data(data, geometry);
%     data = cleanup_neutronics_input_data(data, geometry);
%     [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
%     % Print difference in solution
%     if print_err_bool
%         data.problem.Name = sprintf('z_quad_%s_k2',now_name);
%         outname = ['outputs/',data.problem.Path,'/','ERROR_',data.problem.Name];
%         ye = ExactSol_QuadSolution(DoF.NodeLocations);
%         df = ye - sol.flux{:}; df(abs(df)>1e-12) = 0;
%         write_output_to_vtk_rev2(outname,data,geometry,DoF,df,'flux');
%     end
% end
% clear geometry;
% % Z-Poly Meshes
% % ------------------------------------------------------------------------------
% gname = 'z_mesh_poly_L1_n9_a0.05';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     data.problem.Name = sprintf('SOLUTION_z_poly_%s_k2',now_name);
%     [data, geometry] = process_input_data(data, geometry);
%     data = cleanup_neutronics_input_data(data, geometry);
%     [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
%     % Print difference in solution
%     if print_err_bool
%         data.problem.Name = sprintf('z_poly_%s_k2',now_name);
%         outname = ['outputs/',data.problem.Path,'/','ERROR_',data.problem.Name];
%         ye = ExactSol_QuadSolution(DoF.NodeLocations);
%         df = ye - sol.flux{:}; df(abs(df)>1e-12) = 0;
%         write_output_to_vtk_rev2(outname,data,geometry,DoF,df,'flux');
%     end
% end
% clear geometry;