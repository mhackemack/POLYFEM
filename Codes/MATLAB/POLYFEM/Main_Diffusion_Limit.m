%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Diffusion Limit Runs
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
function Main_Diffusion_Limit()
% Specify some parameters
out_dir = '';
linear_BFs = {'LAGRANGE','MAXENT','PWLD','WACHSPRESS','MV'};
quadratic_BFs = {'LAGRANGE','MAXENT'};
geom_types = {'quad','tri','smooth_poly','rand_poly'};
ep_log_vals = [-1,-2,-3,-4,-5];
% Get Globals, Set Path, and Initialize Domain Space
global glob
glob = get_globals('Office');
fpath = get_path();
addpath(fpath);
addpath([glob.input_path,'Transport']);
[data, ~] = load_user_input();
% Loop through all geometry types
for g=1:length(geom_types)
    geometry = get_geometry(geom_types{g});
    % Run Linear Examples
    for i=1:length(linear_BFs)
        tBF = linear_BFs{i};
        for j=1:length(ep_log_vals)
            data = modify_data(data, 10^(ep_log_vals(j)), tBF);
            [data, geometry] = process_input_data(data, geometry);
            data = cleanup_neutronics_input_data(data, geometry);
            [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
            print_output(out_dir,mesh,DoF,FE,sol.flux,ep_log_vals(j));
        end
    end
    % Run Quadratic Examples
    for i=1:length(quadratic_BFs)
        tBF = quadratic_BFs{i};
        for j=1:length(ep_log_vals)
            data = modify_data(data, 10^(ep_log_vals(j)), tBF);
            [data, geometry] = process_input_data(data, geometry);
            data = cleanup_neutronics_input_data(data, geometry);
            [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
            print_output(out_dir,mesh,DoF,FE,sol.flux,ep_log_vals(j));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = modify_data(data, ep, BF)
data.Neutronics.SpatialMethod = BF;
data.Neutronics.Transport.ScatteringXS = zeros(1,1,1,1);
data.Neutronics.Transport.TotalXS = 1/ep;
data.Neutronics.Transport.AbsorbXS = ep;
data.Neutronics.Transport.ScatteringXS(1,:,:,:) = 1/ep-ep;
data.Neutronics.Transport.ExtSource = ep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = get_geometry(gtype)
global glob
if strcmpi(gtype, 'quad')
    geometry = CartesianGeometry(2,linspace(0,1,21),linspace(0,1,21));
elseif strcmpi(gtype, 'tri')
    tx = linspace(0,1,21);
    [x,y]=meshgrid(tx,tx);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
elseif strcmpi(gtype, 'smooth_poly')
    gname = 'smooth_poly_mesh_L1_n8_a0.15';
    load(strcat(glob.geom_path,gname,'.mat'));
elseif strcmpi(gtype, 'rand_poly')
    gname = 'random_poly_mesh_L1_n4_a0.9';
    load(strcat(glob.geom_path,gname,'.mat'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_output(out_dir,mesh,DoF,FE,flux,ep_val)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%