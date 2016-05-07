%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Transport Linear Solution Run Script
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
% Prepare Project Space
% ------------------------------------------------------------------------------
clc; close all; format long e; clear;
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Office');
inp = 'Transport_Lin_Sol';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, ~] = load_user_input();
% Begin User Input Section
% ------------------------------------------------------------------------------
% BF_names = {'PWLD','WACHSPRESS','MV'};
BF_names = {'WACHSPRESS'};
BF_orders = [1];

data.problem.saveVTKSolution = 1;
data.problem.Dimension = 2;
% End User Input Section
% ------------------------------------------------------------------------------

% Cartesian Meshes
% ------------------------------------------------------------------------------
% x = linspace(0,1,11);
% geometry = CartesianGeometry(2,x,x);
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     for o=1:length(BF_orders)
%         data.Neutronics.FEMDegree = BF_orders(o);
%         if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
%         data.problem.Name = sprintf('cart_%s_k%d',now_name,BF_orders(o));
%         [data, geometry] = process_input_data(data, geometry);
%         data = cleanup_neutronics_input_data(data, geometry);
%         [~, ~, ~, ~, ~] = execute_problem(data, geometry);
%     end
% end
% Triangle Meshes
% ------------------------------------------------------------------------------
% x = linspace(0,1,11); [x,y]=meshgrid(x,x);
% x=x(:); y=y(:);
% tri = delaunayTriangulation(x,y);
% geometry = GeneralGeometry(2, 'Delaunay', tri);
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     for o=1:length(BF_orders)
%         data.Neutronics.FEMDegree = BF_orders(o);
%         if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
%         data.problem.Name = sprintf('tri_%s_k%d',now_name,BF_orders(o));
%         [data, geometry] = process_input_data(data, geometry);
%         data = cleanup_neutronics_input_data(data, geometry);
%         [~, ~, ~, ~, ~] = execute_problem(data, geometry);
%     end
% end
% Shestakov Quad Meshes
% ------------------------------------------------------------------------------
gname = 'shestakov_quad_L1_nc4_a0.25';
load(strcat(glob.geom_path,gname,'.mat'));
for b=1:length(BF_names)
    now_name = BF_names{b};
    data.Neutronics.SpatialMethod = BF_names{b};
    for o=1:length(BF_orders)
        data.Neutronics.FEMDegree = BF_orders(o);
        if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
        data.problem.Name = sprintf('shes_quad_%s_k%d',now_name,BF_orders(o));
        [data, geometry] = process_input_data(data, geometry);
        data = cleanup_neutronics_input_data(data, geometry);
        [~, sol, geometry, DoF, FE] = execute_problem(data, geometry);
    end
end
% Shestakov Poly Meshes
% ------------------------------------------------------------------------------
% gname = 'shestakov_poly_mesh_L1_nc3_a0.15';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     for o=1:length(BF_orders)
%         data.Neutronics.FEMDegree = BF_orders(o);
%         if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
%         data.problem.Name = sprintf('shes_poly_%s_k%d',now_name,BF_orders(o));
%         [data, geometry] = process_input_data(data, geometry);
%         data = cleanup_neutronics_input_data(data, geometry);
%         [~, ~, ~, ~, ~] = execute_problem(data, geometry);
%     end
% end
% Smooth Poly Meshes
% ------------------------------------------------------------------------------
% gname = 'smooth_poly_mesh_L1_n8_a0.15';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     for o=1:length(BF_orders)
%         data.Neutronics.FEMDegree = BF_orders(o);
%         if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
%         data.problem.Name = sprintf('smooth_poly_%s_k%d',now_name,BF_orders(o));
%         [data, geometry] = process_input_data(data, geometry);
%         data = cleanup_neutronics_input_data(data, geometry);
%         [~, ~, ~, ~, ~] = execute_problem(data, geometry);
%     end
% end
% Z-Quad Meshes
% ------------------------------------------------------------------------------
% gname = 'z_mesh_quad_L1_n20_a0.2';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     for o=1:length(BF_orders)
%         data.Neutronics.FEMDegree = BF_orders(o);
%         if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
%         data.problem.Name = sprintf('z_quad_%s_k%d',now_name,BF_orders(o));
%         [data, geometry] = process_input_data(data, geometry);
%         data = cleanup_neutronics_input_data(data, geometry);
%         [~, ~, ~, ~, ~] = execute_problem(data, geometry);
%     end
% end
% Z-Poly Meshes
% ------------------------------------------------------------------------------
% gname = 'z_mesh_poly_L1_n9_a0.05';
% load(strcat(glob.geom_path,gname,'.mat'));
% for b=1:length(BF_names)
%     now_name = BF_names{b};
%     data.Neutronics.SpatialMethod = BF_names{b};
%     for o=1:length(BF_orders)
%         data.Neutronics.FEMDegree = BF_orders(o);
%         if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
%         data.problem.Name = sprintf('z_poly_%s_k%d',now_name,BF_orders(o));
%         [data, geometry] = process_input_data(data, geometry);
%         data = cleanup_neutronics_input_data(data, geometry);
%         [~, ~, ~, ~, ~] = execute_problem(data, geometry);
%     end
% end
% Square Poly Meshes
% ------------------------------------------------------------------------------
gname = 'shestakov_quad_L1_nc4_a0.25';
load(strcat(glob.geom_path,gname,'.mat'));
for b=1:length(BF_names)
    now_name = BF_names{b};
    data.Neutronics.SpatialMethod = BF_names{b};
    for o=1:length(BF_orders)
        data.Neutronics.FEMDegree = BF_orders(o);
        if BF_orders(o) > 1 && ~strcmpi(now_name,'MAXENT'), continue; end
        data.problem.Name = sprintf('shes_quad_%s_k%d',now_name,BF_orders(o));
        [data, geometry] = process_input_data(data, geometry);
        data = cleanup_neutronics_input_data(data, geometry);
        [~, sol, geometry, DoF, FE] = execute_problem(data, geometry);
    end
end