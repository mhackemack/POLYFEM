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
% ---
geom_in.Dimension = 2;
geom_in.GeometryType = 'tri';
pnum = [4,16,64,256,1024,4096,16384,65536];
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
sdm = {'PWLD','WACHSPRESS','MV','MAXENT'};
fedeg = [1,2];
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
dat_in.refineMesh = 1;
dat_in.refinementLevels = 7;
dat_in.AMRIrregularity = 1;
dat_in.refinementTolerance = 0.0;
dat_in.projectSolution = 0;
% ---
sigt = [1,10,100,1000];
dat_in.TotalXS = 1;
dat_in.RHSFunc = {@ZeroTransportFunction};
dat_in.SolFunc = {@ExactSol_IncidentLeftFace_2D_45degDown_LS4_sigt10};
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
% [data, geometry] = load_user_input(dat_in, geom_in);
% Modify path
if strcmpi(geom_in.GeometryType, 'cart')
    gt = 'Cartesian';
elseif strcmpi(geom_in.GeometryType, 'tri')
    gt = 'Triangular';
elseif strcmpi(geom_in.GeometryType, 'poly')
    gt = 'Polygonal';
end
problem_path = ['Transport/PureAbsorber/',prob_name,'/',gt];
% data.problem.Name = dat_in.Name;
% Execute Problem Suite
% ------------------------------------------------------------------------------
% Run polygonal problem
if strcmpi(geom_in.GeometryType, 'poly')
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
                % Run problem iteration
                msg = sprintf('fedeg: %d of %d, sdm: %d of %d, sigt: %d of %d.',k,length(fedeg),s,length(sdm),t,length(sigt));
                disp(msg);
                % Loop through poly geometries
                for g=1:length(pnum)
                    geom_in.PolyNum = pnum(g);
                    % Load data and build problem name
                    [data, geometry] = load_user_input(dat_in, geom_in);
                    [data, geometry] = process_input_data(data, geometry);
                    data = cleanup_neutronics_input_data(data, geometry);
                    [data, sol, geometry, DoF, ~] = execute_problem(data, geometry);
                    cell_vols = (geometry.CellVolume).^(1/geom_in.Dimension);
                    vol_ave = mean(cell_vols);
                    vol_max = max(cell_vols);
                end
            end
        end
    end
% Run triangular or cartesian problem
else
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
                % Load data and build problem name
                [data, geometry] = load_user_input(dat_in, geom_in);
                data.problem.Name = sprintf('%s%d_sigt%d',sdm{s},fedeg(k),sigt(t));
                data.problem.Path = problem_path;
                % Run problem iteration
                msg = sprintf('fedeg: %d of %d, sdm: %d of %d, sigt: %d of %d.',k,length(fedeg),s,length(sdm),t,length(sigt));
                disp(msg);
                [data, geometry] = process_input_data(data, geometry);
                data = cleanup_neutronics_input_data(data, geometry);
                [data, sol, geometry, ~, ~] = execute_problem(data, geometry);
                [dofs,err] = retrieve_MMS_error(data,sol);
                % Build some output structures
                totcells = zeros(dat_in.refinementLevels+1,1);
                avecellvol = zeros(dat_in.refinementLevels+1,1);
                maxcellvol = zeros(dat_in.refinementLevels+1,1);
                totflux    = zeros(dat_in.refinementLevels+1,1);
                % Collect results and output
                for r=1:(dat_in.refinementLevels+1)
                    rr = r-1;
                    totcells(r)   = sum(sol{r}.CellVertexNumbers);
                    avecellvol(r) = sol{r}.AverageCellMeasure;
                    maxcellvol(r) = sol{r}.MaxCellMeasure;
                    totflux(r)    = sol{r}.TotalMaterialFlux;
                end
                % Output the results
                dlmwrite(['outputs/',data.problem.Path,'/',data.problem.Name,'_outdata.dat'],[totcells,dofs,avecellvol,maxcellvol,totflux,err]);
            end
        end
    end
end
