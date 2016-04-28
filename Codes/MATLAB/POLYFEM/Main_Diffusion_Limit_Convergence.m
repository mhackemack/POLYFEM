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
function Main_Diffusion_Limit_Convergence(out_dir)
clearvars -except out_dir; close all; clc; format long e;
% Specify some parameters
linear_BFs = {'PWLD'};
quadratic_BFs = {'MAXENT'};
geom_types = {'quad'};
% linear_BFs = {'WACHSPRESS','PWLD','MV','MAXENT'};
% quadratic_BFs = {'WACHSPRESS','PWLD','MV','MAXENT'};
% geom_types = {'quad','Sq_poly'};
% ep_log_vals = [-6];
ep_log_vals = [0,-1,-2,-3];
% ep_log_vals = [-1,-2,-3,-4,-5,-6];
% Get Globals, Set Path, and Initialize Domain Space
global glob
glob = get_globals('Office');
fpath = get_path();
addpath(fpath);
addpath([glob.input_path,'DL_Transport']);
transdata = load_user_input();
addpath([glob.input_path,'DL_Diffusion']);
diffdata = load_user_input();
% Allocate memory space
lin_sol_err  = zeros(length(ep_log_vals), length(linear_BFs), length(geom_types));
quad_sol_err = zeros(length(ep_log_vals), length(quadratic_BFs), length(geom_types));
% Loop through all geometry types
for g=1:length(geom_types)
    geometry = get_geometry(geom_types{g});
    % Run Linear Cases
    tdata = transdata; tdata.Neutronics.FEMDegree = 1;
    ddata = diffdata; ddata.Neutronics.FEMDegree = 1;
    % Loop through basis functions
    for i=1:length(linear_BFs)
        tBF = linear_BFs{i};
        if ~check_BF_geom_combo(tBF, geometry), continue; end
        % Diffusion problem
        ddata = modify_diffusion_data(ddata, 1, tBF);
        [ddata, geometry] = process_input_data(ddata, geometry);
        ddata = cleanup_neutronics_input_data(ddata, geometry);
        [ddata, dsol, ~, DoF, FE] = execute_problem(ddata, geometry);
        % Transport problems
        for j=1:length(ep_log_vals)
            tdata = modify_transport_data(tdata, 10^(ep_log_vals(j)), tBF);
            [tdata, geometry] = process_input_data(tdata, geometry);
            tdata = cleanup_neutronics_input_data(tdata, geometry);
            [tdata, tsol, ~, ~, ~] = execute_problem(tdata, geometry);
            lin_sol_err(j,i,g) = calc_diff_trans_error(geometry,DoF,FE,dsol.flux{1},tsol.flux{1});
        end
    end
    % Run Quadratic Cases
    tdata = transdata; tdata.Neutronics.FEMDegree = 2;
    ddata = diffdata; ddata.Neutronics.FEMDegree = 2;
    % Loop through basis functions
    for i=1:length(quadratic_BFs)
        tBF = quadratic_BFs{i};
        if ~check_BF_geom_combo(tBF, geometry), continue; end
        % Diffusion problem
        ddata = modify_diffusion_data(ddata, 1, tBF);
        [ddata, geometry] = process_input_data(ddata, geometry);
        ddata = cleanup_neutronics_input_data(ddata, geometry);
        [ddata, dsol, ~, DoF, FE] = execute_problem(ddata, geometry);
        % Transport problems
        for j=1:length(ep_log_vals)
            tdata = modify_transport_data(tdata, 10^(ep_log_vals(j)), tBF);
            [tdata, geometry] = process_input_data(tdata, geometry);
            tdata = cleanup_neutronics_input_data(tdata, geometry);
            [tdata, tsol, ~, ~, ~] = execute_problem(tdata, geometry);
            quad_sol_err(j,i,g) = calc_diff_trans_error(geometry,DoF,FE,dsol.flux{1},tsol.flux{1});
        end
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = modify_transport_data(data, ep, BF)
data.Neutronics.SpatialMethod = BF;
data.Neutronics.Transport.ScatteringXS = zeros(1,1,1,1);
data.Neutronics.Transport.TotalXS = 1/ep;
data.Neutronics.Transport.AbsorbXS = ep;
data.Neutronics.Transport.ScatteringXS(1,:,:,:) = 1/ep-ep;
data.Neutronics.Transport.ExtSource = ep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = modify_diffusion_data(data, ep, BF)
data.Neutronics.SpatialMethod = BF;
% data.Neutronics.Diffusion.ScatteringXS = zeros(1,1,1,1);
% data.Neutronics.Diffusion.TotalXS = ep;
% data.Neutronics.Diffusion.DiffXS = ep/3;
% data.Neutronics.Diffusion.AbsorbXS = 0;
% data.Neutronics.Diffusion.ScatteringXS(1,:,:,:) = 0;
% data.Neutronics.Diffusion.ExtSource = ep;
data.Neutronics.Diffusion.ScatteringXS = zeros(1,1,1,1);
data.Neutronics.Diffusion.TotalXS = 1;
data.Neutronics.Diffusion.DiffXS = 1/3;
data.Neutronics.Diffusion.AbsorbXS = 0;
data.Neutronics.Diffusion.ScatteringXS(1,:,:,:) = 0;
data.Neutronics.Diffusion.ExtSource = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = get_geometry(gtype)
global glob
if strcmpi(gtype, 'quad')
%     geometry = CartesianGeometry(2,linspace(0,1,11),linspace(0,1,11));
    geometry = CartesianGeometry(2,linspace(0,1,21),linspace(0,1,21));
elseif strcmpi(gtype, 'tri')
    tx = linspace(0,1,11);
    [x,y]=meshgrid(tx,tx);
    x=x(:);y=y(:);
    tri = delaunayTriangulation(x,y);
    geometry = GeneralGeometry(2, 'Delaunay', tri);
elseif strcmpi(gtype, 'smooth_poly')
    gname = 'smooth_poly_mesh_L1_n8_a0.15';
    load(strcat(glob.geom_path,gname,'.mat'));
elseif strcmpi(gtype, 'Sq_poly')
    gname = 'PolyMesh_SqDomain_L1_n256';
    load(strcat(glob.geom_path,gname,'.mat'));
elseif strcmpi(gtype, 'z-poly')
    gname = 'z_mesh_poly_L1_n9_a0.05';
    load(strcat(glob.geom_path,gname,'.mat'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = check_BF_geom_combo(tBF, geometry)
gtype = geometry.MeshType;
if strcmpi(gtype, 'Quadrilateral') || strcmpi(gtype, 'Triangle')
    out = true;
else
    if strcmpi(tBF,'lagrange') || strcmpi(tBF,'serendipity')
        out = false;
    else
        out = true;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = calc_diff_trans_error(mesh,DoF,FE,dsol,tsol)
out = 0;
% Loop through cells
for c=1:mesh.TotalCells
    cf = mesh.CellFaces{c};
    cfflags = mesh.FaceID(cf); cfflags = cfflags(cfflags~=0);
    if isempty(cfflags)
        cdofs = DoF.ConnectivityArray{c};
        M = FE.CellMassMatrix{c};
        dval = dsol(cdofs);
        tval = tsol(cdofs);
        out = out + (M*(dval-tval))'*(dval-tval);
    end
end
out = sqrt(out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%