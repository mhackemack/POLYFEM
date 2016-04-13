%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D Pure Absorber Problem Run Script
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
global glob; glob = get_globals('Home');
inp = 'Transport_PureAbsorber';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
prob_name = 'IncidentLeftFace_2D_45degDown_LS4';
cart_run_bool = true;
tri_run_bool  = true;
poly_run_bool = false;
split_poly_run_bool = false;
% ---
geom_in.Dimension = 2;
geom_in.GeometryType = 'tri';
% pnum = [16,64,256,1024,4096,16384];
pnum = [262144];
% geom_in.PolyNum = [4,16,64,256,1024,4096,16384,65536];
geom_in.Lx = 1; geom_in.ncellx = 4;
geom_in.Ly = 1; geom_in.ncelly = 4;
geom_in.Lz = 1; geom_in.ncellz = 4;
geom_in.xmin_bound_type = glob.Function;
geom_in.xmax_bound_type = glob.Vacuum;
geom_in.ymin_bound_type = glob.Vacuum;
geom_in.ymax_bound_type = glob.Vacuum;
geom_in.zmin_bound_type = glob.Reflecting;
geom_in.zmax_bound_type = glob.Reflecting;
geom_in.xmin_val = @BoundaryFunc_IncidentLeftFace_2D_45degDown_LS4;
geom_in.xmax_val = 0;
geom_in.ymin_val = 0;
geom_in.ymax_val = 0;
geom_in.zmin_val = 0;
geom_in.zmax_val = 0;
% ---
% sdm = {'PWLD'};
sdm = {'PWLD','WACHSPRESS','MV','MAXENT'};
fedeg = [1];
dat_in.SpatialMethod = 'PWLD';
dat_in.FEMDegree = 2;
dat_in.FEMLumping = false;
% ---
dat_in.QuadType = 'LS';
dat_in.SnLevels = 4;
dat_in.AzimuthalLevels = 24;
dat_in.PolarLevels = 4;
dat_in.QuadAngles  = [1,-1]/norm([1,-1]);  % Angles for manual set
dat_in.QuadWeights = 1;                  % Weights for manual set
% ---
dat_in.refineMesh = 0;
dat_in.refinementLevels = 0;
dat_in.AMRIrregularity = 1;
dat_in.refinementTolerance = 0.0;
dat_in.projectSolution = 0;
% ---
sigt = 50;
% sigt = [1,10,100,1000];
dat_in.TotalXS = 1;
dat_in.RHSFunc = {@ZeroTransportFunction};
dat_in.SolFunc = {@ExactSol_IncidentLeftFace_2D_45degDown_LS4_sigt10};
% Execute Problem Suite
% ------------------------------------------------------------------------------
print_heading(now, date);
pmax = log2(sqrt(max(pnum)));
% Run cartesian problem
% ------------------------------------------------------------------------------
if cart_run_bool
    geom_in.GeometryType = 'cart';
    problem_path = ['Transport/PureAbsorber/',prob_name,'/Cartesian'];
    % Loop through FEM degrees
    for k=1:length(fedeg)
        dat_in.FEMDegree = fedeg(k);
        % Loop through basis functions
        for s=1:length(sdm)
            dat_in.SpatialMethod = sdm{s};
            % Loop through total cross sections
            for t=1:length(sigt)
                dat_in.TotalXS = sigt(t);
                dat_in.SolFunc{1} = str2func(['ExactSol_',prob_name,'_sigt',num2str(sigt(t))]);
                % Build some output structures
                totcells = zeros(length(pnum),1);
                avecellvol = zeros(length(pnum),1);
                maxcellvol = zeros(length(pnum),1);
                totflux    = zeros(length(pnum),1);
                dofs       = zeros(length(pnum),1);
                err        = zeros(length(pnum),1);
                % Loop through meshes
                for g=1:length(pnum)
                    msg = sprintf('fedeg: %d of %d, sdm: %d of %d, sigt: %d of %d, n: %d of %d.',k,length(fedeg),s,length(sdm),t,length(sigt),g,length(pnum));
                    disp(msg);
                    geom_in.ncellx = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncelly = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncellz = (pnum(g))^(1/geom_in.Dimension);
                    % Load data and build problem name
                    data = load_user_input(dat_in, geom_in);
                    [data,geometry] = load_geometry_input(data, geom_in);
                    data.problem.Name = sprintf('%s%d_sigt%d_n%d',sdm{s},fedeg(k),sigt(t),pnum(g));
                    data.problem.Path = problem_path;
                    % Run problem iteration
                    [data, geometry] = process_input_data(data, geometry);
                    data = cleanup_neutronics_input_data(data, geometry);
                    [data, sol, ~, ~, ~] = execute_problem(data, geometry);
                    [tdofs,terr] = retrieve_MMS_error(data,sol);
                    % Collect results and output
                    totcells(g)   = sum(sol.CellVertexNumbers);
                    avecellvol(g) = sol.AverageCellMeasure;
                    maxcellvol(g) = sol.MaxCellMeasure;
                    totflux(g)    = sol.TotalMaterialFlux;
                    dofs(g)       = tdofs;
                    err(g)        = terr;
                end
                % Output the results
                data.problem.Name = sprintf('%s%d_sigt%d',sdm{s},fedeg(k),sigt(t));
                dlmwrite(['outputs/',data.problem.Path,'/',data.problem.Name,'_n',num2str(pmax),'_outdata.dat'],[totcells,dofs,avecellvol,maxcellvol,totflux,err]);
            end
        end
    end
end
% Run triangular problem
% ------------------------------------------------------------------------------
if tri_run_bool
    geom_in.GeometryType = 'tri';
    problem_path = ['Transport/PureAbsorber/',prob_name,'/Triangular'];
    % Build triangular geometries and store them - the construction can get long
    % for the higher cell count cases.
    tri_geoms = cell(length(pnum),1);
    data = load_user_input(dat_in, geom_in);
    for g=1:length(pnum)
        geom_in.ncellx = (pnum(g))^(1/geom_in.Dimension);
        geom_in.ncelly = (pnum(g))^(1/geom_in.Dimension);
        geom_in.ncellz = (pnum(g))^(1/geom_in.Dimension);
        [data,tri_geoms{g}] = load_geometry_input(data, geom_in);
    end
    % Loop through FEM degrees
    for k=1:length(fedeg)
        dat_in.FEMDegree = fedeg(k);
        % Loop through basis functions
        for s=1:length(sdm)
            dat_in.SpatialMethod = sdm{s};
            % Loop through total cross sections
            for t=1:length(sigt)
                dat_in.TotalXS = sigt(t);
                dat_in.SolFunc{1} = str2func(['ExactSol_',prob_name,'_sigt',num2str(sigt(t))]);
                % Build some output structures
                totcells = zeros(length(pnum),1);
                avecellvol = zeros(length(pnum),1);
                maxcellvol = zeros(length(pnum),1);
                totflux    = zeros(length(pnum),1);
                dofs       = zeros(length(pnum),1);
                err        = zeros(length(pnum),1);
                % Loop through meshes
                for g=1:length(pnum)
                    msg = sprintf('fedeg: %d of %d, sdm: %d of %d, sigt: %d of %d, n: %d of %d.',k,length(fedeg),s,length(sdm),t,length(sigt),g,length(pnum));
                    disp(msg);
                    geom_in.ncellx = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncelly = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncellz = (pnum(g))^(1/geom_in.Dimension);
                    % Load data and build problem name
                    data = load_user_input(dat_in, geom_in);
                    data = set_PA_BC_flags(data, geom_in);
                    data.problem.Name = sprintf('%s%d_sigt%d_n%d',sdm{s},fedeg(k),sigt(t),pnum(g));
                    data.problem.Path = problem_path;
                    % Run problem iteration
                    [data, ~] = process_input_data(data, tri_geoms{g});
                    data = cleanup_neutronics_input_data(data, tri_geoms{g});
                    [data, sol, ~, ~, ~] = execute_problem(data, tri_geoms{g});
                    [tdofs,terr] = retrieve_MMS_error(data,sol);
                    % Collect results and output
                    totcells(g)   = sum(sol.CellVertexNumbers);
                    avecellvol(g) = sol.AverageCellMeasure;
                    maxcellvol(g) = sol.MaxCellMeasure;
                    totflux(g)    = sol.TotalMaterialFlux;
                    dofs(g)       = tdofs;
                    err(g)        = terr;
                end
                % Output the results
                data.problem.Name = sprintf('%s%d_sigt%d',sdm{s},fedeg(k),sigt(t));
                dlmwrite(['outputs/',data.problem.Path,'/',data.problem.Name,'_n',num2str(pmax),'_outdata.dat'],[totcells,dofs,avecellvol,maxcellvol,totflux,err]);
            end
        end
    end
    clear tri_geoms;
end
% Run polygonal problem
% ------------------------------------------------------------------------------
if poly_run_bool
    geom_in.GeometryType = 'poly';
    problem_path = ['Transport/PureAbsorber/',prob_name,'/Polygonal'];
    % Loop through FEM degrees
    for k=1:length(fedeg)
        dat_in.FEMDegree = fedeg(k);
        % Loop through basis functions
        for s=1:length(sdm)
            dat_in.SpatialMethod = sdm{s};
            % Loop through total cross sections
            for t=1:length(sigt)
                dat_in.TotalXS = sigt(t);
                dat_in.SolFunc{1} = str2func(['ExactSol_',prob_name,'_sigt',num2str(sigt(t))]);
                % Build some output structures
                totcells = zeros(length(pnum),1);
                avecellvol = zeros(length(pnum),1);
                maxcellvol = zeros(length(pnum),1);
                totflux    = zeros(length(pnum),1);
                dofs       = zeros(length(pnum),1);
                err        = zeros(length(pnum),1);
                % Loop through meshes
                for g=1:length(pnum)
                    msg = sprintf('fedeg: %d of %d, sdm: %d of %d, sigt: %d of %d, n: %d of %d.',k,length(fedeg),s,length(sdm),t,length(sigt),g,length(pnum));
                    disp(msg);
                    geom_in.PolyNum = pnum(g);
                    geom_in.ncellx = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncelly = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncellz = (pnum(g))^(1/geom_in.Dimension);
                    % Load data and build problem name
                    data = load_user_input(dat_in, geom_in);
                    [data,geometry] = load_geometry_input(data, geom_in);
                    data.problem.Name = sprintf('%s%d_sigt%d_n%d',sdm{s},fedeg(k),sigt(t),pnum(g));
                    data.problem.Path = problem_path;
                    % Run problem iteration
                    [data, geometry] = process_input_data(data, geometry);
                    data = cleanup_neutronics_input_data(data, geometry);
                    [data, sol, ~, ~, ~] = execute_problem(data, geometry);
                    [tdofs,terr] = retrieve_MMS_error(data,sol);
                    % Collect results and output
                    totcells(g)   = sum(sol.CellVertexNumbers);
                    avecellvol(g) = sol.AverageCellMeasure;
                    maxcellvol(g) = sol.MaxCellMeasure;
                    totflux(g)    = sol.TotalMaterialFlux;
                    dofs(g)       = tdofs;
                    err(g)        = terr;
                end
                % Output the results
                data.problem.Name = sprintf('%s%d_sigt%d',sdm{s},fedeg(k),sigt(t));
                dlmwrite(['outputs/',data.problem.Path,'/',data.problem.Name,'_n',num2str(pmax),'_outdata.dat'],[totcells,dofs,avecellvol,maxcellvol,totflux,err]);
            end
        end
    end
end
% Run split-polygonal problem
% ------------------------------------------------------------------------------
if split_poly_run_bool
    geom_in.GeometryType = 'split_poly';
    problem_path = ['Transport/PureAbsorber/',prob_name,'/SplitPolygonal'];
    % Loop through FEM degrees
    for k=1:length(fedeg)
        dat_in.FEMDegree = fedeg(k);
        % Loop through basis functions
        for s=1:length(sdm)
            dat_in.SpatialMethod = sdm{s};
            % Loop through total cross sections
            for t=1:length(sigt)
                dat_in.TotalXS = sigt(t);
                dat_in.SolFunc{1} = str2func(['ExactSol_',prob_name,'_sigt',num2str(sigt(t))]);
                % Build some output structures
                totcells = zeros(length(pnum),1);
                avecellvol = zeros(length(pnum),1);
                maxcellvol = zeros(length(pnum),1);
                totflux    = zeros(length(pnum),1);
                dofs       = zeros(length(pnum),1);
                err        = zeros(length(pnum),1);
                % Loop through meshes
                for g=1:length(pnum)
                    msg = sprintf('fedeg: %d of %d, sdm: %d of %d, sigt: %d of %d, n: %d of %d.',k,length(fedeg),s,length(sdm),t,length(sigt),g,length(pnum));
                    disp(msg);
                    geom_in.PolyNum = pnum(g);
                    geom_in.ncellx = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncelly = (pnum(g))^(1/geom_in.Dimension);
                    geom_in.ncellz = (pnum(g))^(1/geom_in.Dimension);
                    % Load data and build problem name
                    data = load_user_input(dat_in, geom_in);
                    [data,geometry] = load_geometry_input(data, geom_in);
                    data.problem.Name = sprintf('%s%d_sigt%d_n%d',sdm{s},fedeg(k),sigt(t),pnum(g));
                    data.problem.Path = problem_path;
                    % Run problem iteration
                    [data, geometry] = process_input_data(data, geometry);
                    data = cleanup_neutronics_input_data(data, geometry);
                    [data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
                    [tdofs,terr] = retrieve_MMS_error(data,sol);
                    % Collect results and output
                    totcells(g)   = sum(sol.CellVertexNumbers);
                    avecellvol(g) = sol.AverageCellMeasure;
                    maxcellvol(g) = sol.MaxCellMeasure;
                    totflux(g)    = sol.TotalMaterialFlux;
                    dofs(g)       = tdofs;
                    err(g)        = terr;
                end
                % Output the results
                data.problem.Name = sprintf('%s%d_sigt%d',sdm{s},fedeg(k),sigt(t));
                dlmwrite(['outputs/',data.problem.Path,'/',data.problem.Name,'_n',num2str(pmax),'_outdata.dat'],[totcells,dofs,avecellvol,maxcellvol,totflux,err]);
            end
        end
    end
end