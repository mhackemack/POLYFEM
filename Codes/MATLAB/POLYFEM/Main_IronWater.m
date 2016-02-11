%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          IronWater Problem Run Script
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
clc; close all; format long e
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Home');
inp = 'IronWater';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
dat_in.GeometryType = 'cart';
dat_in.SpatialMethod = 'MAXENT';
dat_in.FEMDegree = 2;
dat_in.FEMLumping = false;
% ---
dat_in.QuadType = 'LS';
dat_in.SnLevels = 4;
dat_in.AzimuthalLevels = 24;
dat_in.PolarLevels = 1;
% ---
dat_in.refinementLevels = 27;
dat_in.AMRIrregularity = 3;
dat_in.refinementTolerance = 1/10;
dat_in.projectSolution = 0;
% ---
dat_in.DSASolveMethod = 'PCG';
dat_in.DSAPreconditioner = 'eisenstat';
dat_in.DSATolerance = 1e-3;
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input(dat_in);
% Modify path
if dat_in.projectSolution == 1
    bs = 'WITHBootstrapping';
elseif dat_in.projectSolution == 0
    bs = 'NOBootstrapping';
end
if strcmpi(dat_in.GeometryType, 'cart')
    gt = 'Cartesian';
elseif strcmpi(dat_in.GeometryType, 'tri')
    gt = 'Triangular';
end
data.problem.Path = ['Transport/IronWater/',bs,'/',gt,'/',dat_in.SpatialMethod,num2str(dat_in.FEMDegree)];
% Modify name
if strcmpi(dat_in.QuadType, 'LS')
    oname = [dat_in.QuadType,num2str(dat_in.SnLevels)];
elseif strcmpi(dat_in.QuadType, 'PGLC')
    oname = [dat_in.QuadType,num2str(dat_in.AzimuthalLevels),'-',num2str(dat_in.PolarLevels)];
end
if abs(dat_in.refinementTolerance) < 1e-8
    oname = [oname, '_uniform'];
else
    oname = [oname,sprintf('_Irr=%d_tol=%1.3f',dat_in.AMRIrregularity,dat_in.refinementTolerance)];
end
oname = [oname,'_',dat_in.DSASolveMethod,'_',dat_in.DSAPreconditioner];
data.problem.Name = oname;
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, sol, ~, ~, ~] = execute_problem(data, geometry);
% Build data storage structures
nr = data.problem.refinementLevels + 1;
dofnum = zeros(nr, 1); maxv = 0; 
cell_verts = cell(nr, 1); cell_nums = zeros(nr, 1);
ave_flux = zeros(nr, 3); tot_flux = zeros(nr, 3);
tot_int = zeros(nr, 3); tot_abs = zeros(nr, 3);
out_iters = zeros(nr, 1); out_DSA_iters = zeros(nr, 1);
out_times = zeros(nr, 1);
% Loop through AMR cycles
ddir = ['outputs/',data.problem.Path,'/',data.problem.Name];
for rlvl=0:nr-1
    (fprintf(1,'Refinement Calculation: %d of %d.\n',rlvl,nr-1));
    cname = ['_',num2str(rlvl)];
    % Load data structures
    load([ddir,'_geometry',cname,'.mat']);
    load([ddir,'_sol',cname,'.mat']);
    dofnum(rlvl+1) = sol.SpatialDoFs;
    cell_nums(rlvl+1) = geometry.TotalCells;
    iter = sol.iter;
    DSA_iters = (sol.DSA_iters); num_DSA_iters = sum(DSA_iters);
    times = sol.times; tot_time = sum(times); 
    out_times(rlvl+1) = tot_time;
    out_iters(rlvl+1) = iter;
    out_DSA_iters(rlvl+1) = num_DSA_iters;
    cell_verts{rlvl+1} = sol.CellVertexNumbers;
    tmaxv = length(sol.CellVertexNumbers);
    if tmaxv > maxv, maxv = tmaxv; end
    ave_flux(rlvl+1,:) = sol.AverageMaterialFlux';
    tot_flux(rlvl+1,:) = sol.TotalMaterialFlux';
    clear geometry;
end
out_cell_verts = zeros(nr, maxv);
for rlvl=0:nr-1
    nv = length(sol.CellVertexNumbers);
    out_cell_verts(rlvl+1,1:nv) = sol.CellVertexNumbers;
end
% Print outputs
% ------------------------------------------------------------------------------
dlmwrite([ddir,'_numcellverts','.dat'],out_cell_verts);
dlmwrite([ddir,'_outdata','.dat'],[cell_nums,dofnum,out_iters,out_DSA_iters,out_times,tot_flux,ave_flux],'precision','%18.14e');

