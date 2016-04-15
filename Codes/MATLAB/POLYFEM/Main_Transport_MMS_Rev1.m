%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Transport MMS Script (Rev1)
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
% Clear Project Space
% ------------------------------------------------------------------------------
clc; close all; clear; format long e;
fpath = get_path(); 
addpath(fpath);
% Populate global space
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Home');
glob.print_info = true;
inp = 'Transport_MMS_Rev1';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Being User Input Section
% ------------------------------------------------------------------------------
path = 'Transport_MMS/Gauss2D';
dat_in.prob = 'gauss_iso';
% ---
geom_in.Dimension = 2;
geom_in.GeometryType = 'cart';
geom_in.Lx = 1; geom_in.ncellx = 4;
geom_in.Ly = 1; geom_in.ncelly = 4;
% ---
fedeg = [1];
sdm = {'PWLD'};
% fedeg = [1,2];
% sdm = {'PWLD','MV','MAXENT'};
% ---
dat_in.lvls = 24;
dat_in.irr = 1;
dat_in.tol = 0.1;
% Execute Problem Suite
% ------------------------------------------------------------------------------
print_heading(now, date);
[data,geometry] = load_user_input(dat_in, geom_in);
% Loop through finite element order
for k=1:length(fedeg)
    data.Neutronics.FEMDegree = fedeg(k);
    % Loop through basis functions
    for s=1:length(sdm)
        data.Neutronics.SpatialMethod = sdm{s};
        data.problem.Path = sprintf('%s/%s_k%d',path,sdm{s},fedeg(k));
        data.problem.Name = sprintf('%s_Irr=%d_tol=%g',geom_in.GeometryType,dat_in.irr,dat_in.tol);
        [data, geometry] = process_input_data(data, geometry);
        data = cleanup_neutronics_input_data(data, geometry);
        [data, sol, ~, ~, ~] = execute_problem(data, geometry);
        [dofs,err] = retrieve_MMS_error(data,sol);
        dlmwrite(['outputs/',data.problem.Path,'/',data.problem.Name,'_mmserror','.dat'],[dofs,err],'precision','%18.14e');
        maxv = 0;
        for i=0:data.problem.refinementLevels
            if length(sol{i+1}.CellVertexNumbers) > maxv
                maxv = length(sol{i+1}.CellVertexNumbers);
            end
        end
        cellverts = zeros(data.problem.refinementLevels+1, maxv);
        for i=0:data.problem.refinementLevels
            nv = length(sol{i+1}.CellVertexNumbers);
            cellverts(i+1,1:nv) = sol{i+1}.CellVertexNumbers;
        end
        dlmwrite(['outputs/',data.problem.Path,'/',data.problem.Name,'_numcellverts','.dat'],cellverts);
    end
end

