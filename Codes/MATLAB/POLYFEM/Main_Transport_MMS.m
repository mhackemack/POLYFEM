%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Transport Method of Manufactured Solutions (MMS) Script
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
clc; close all; format long e
fpath = get_path(); 
addpath(fpath);
% Populate global space
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Home');
glob.print_info = true;
inp = 'Transport_MMS';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Load generic data
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input();
data.problem.Path = sprintf('%s/%s_k%d',data.problem.Path,data.Neutronics.SpatialMethod,data.Neutronics.FEMDegree);
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
[data, sol, geometry, DoF, FE] = execute_problem(data, geometry);
% Process Outputs
% ------------------------------------------------------------------------------
ddir = ['outputs/',data.problem.Path,'/',data.problem.Name];
[dofs,err] = retrieve_MMS_error(data,sol);
dlmwrite([ddir,'_mmserror','.dat'],[dofs,err],'precision','%18.14e');
% Loop through refinements to get max max number of vertices
if ~data.problem.refineMesh || data.problem.refinementLevels == 0
    dlmwrite([ddir,'_numcellverts','.dat'],sol.CellVertexNumbers);
else
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
    dlmwrite([ddir,'_numcellverts','.dat'],cellverts);
end
