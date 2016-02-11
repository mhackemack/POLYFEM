
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Searchlight Problem Run Script
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
inp = 'Searchlight';
addpath([glob.input_path,inp]); % This one must be last to properly switch input files
% Load data and perform error checking
% ------------------------------------------------------------------------------
print_heading(now, date);
[data, geometry] = load_user_input();
% Modify path
data.problem.Path = sprintf('Transport/Searchlight/%s_k%d',data.Neutronics.SpatialMethod,data.Neutronics.FEMDegree);
[data, geometry] = process_input_data(data, geometry);
data = cleanup_neutronics_input_data(data, geometry);
% Execute Problem Suite
% ------------------------------------------------------------------------------
[data, ~, ~, ~, ~] = execute_problem(data, geometry);
% Process Outputs
% ------------------------------------------------------------------------------
adir = data.Neutronics.Transport.QuadAngles';
% Build data storage structures
nr = data.problem.refinementLevels + 1;
dofnum = zeros(nr,1);
lfluxsum = zeros(nr,1); rfluxsum = zeros(nr,1);
lfluxbndsum = zeros(nr,1); rfluxbndsum = zeros(nr,1);
influx = data.Neutronics.Transport.BCVals{2}*dot([1,0], data.Neutronics.Transport.QuadAngles)*0.2;
lfc = cell(nr,1); rfc = cell(nr,1);
lpos = cell(nr,1); rpos = cell(nr,1);
lflux = cell(nr,1); rflux = cell(nr,1);
% Loop through AMR cycles
ddir = ['outputs/',data.problem.Path,'/',data.problem.Name];
for rlvl=0:nr-1
    (fprintf(1,'Refinement Calculation: %d of %d.\n',rlvl,nr-1));
    cname = ['_',num2str(rlvl)];
    % Load data structures
    load([ddir,'_data',cname,'.mat']);
    load([ddir,'_geometry',cname,'.mat']);
    load([ddir,'_DoF',cname,'.mat']);
    load([ddir,'_FE',cname,'.mat']);
    load([ddir,'_sol',cname,'.mat']);
    flux = sol.flux{:};
    dofnum(rlvl+1) = DoF.TotalDoFs;
    lcounter = 1; rcounter = 1;
    tlpos = []; trpos = [];
    tlflux = []; trflux = [];
    % Loop through boundary faces
    for ff=1:geometry.TotalBoundaryFaces
        f = geometry.BoundaryFaces(ff);
        fdir = geometry.FaceNormal(f,:);
        fnorm = fdir*adir;
        fc = geometry.FaceCenter(f,:);
        fcn = DoF.FaceCellNodes{f,1};
        % Appropriate left boundary
        if abs(fc(1)) < 1e-14
            pos = DoF.NodeLocations(fcn,2); af = flux(fcn);
            M = FE.FaceMassMatrix{f};
            lfluxsum(rlvl+1) = lfluxsum(rlvl+1) + sum(fnorm*M*flux(fcn));
            lfc{rlvl+1} = [lfc{rlvl+1}; fc(2)];
            [pos,ind] = sort(pos);
            tlpos{lcounter} = pos; tlflux{lcounter} = af(ind);
            lcounter = lcounter + 1;
            % Within entrance bounds
            if fc(2) > 0.2 && fc(2) < 0.4
                lfluxbndsum(rlvl+1) = lfluxbndsum(rlvl+1) + sum(fnorm*M*flux(fcn));
            end
        end
        % Appropriate right boundary
        if abs(fc(1)-1.0) < 1e-14
            pos = DoF.NodeLocations(fcn,2); af = flux(fcn);
            M = FE.FaceMassMatrix{f};
            rfluxsum(rlvl+1) = rfluxsum(rlvl+1) + sum(fnorm*M*flux(fcn));
            rfc{rlvl+1} = [rfc{rlvl+1}; fc(2)];
            [pos,ind] = sort(pos);
            trpos{rcounter} = pos; trflux{rcounter} = af(ind);
            rcounter = rcounter + 1;
            % Within exit bounds
            if fc(2) > 0.6 && fc(2) < 0.8
                rfluxbndsum(rlvl+1) = rfluxbndsum(rlvl+1) + sum(fnorm*M*flux(fcn));
            end
        end
    end
    % Organize left boundary outputs
    [~,ind] = sort(lfc{rlvl+1});
    for j=1:length(ind)
        lpos{rlvl+1} = [lpos{rlvl+1};tlpos{ind(j)}];
        lflux{rlvl+1} = [lflux{rlvl+1};tlflux{ind(j)}];
    end
    dlmwrite([ddir,'_leftflux',cname,'.dat'],[lpos{rlvl+1},lflux{rlvl+1}],'precision','%18.14e');
    % Organize right boundary outputs
    [~,ind] = sort(rfc{rlvl+1});
    for j=1:length(ind)
        rpos{rlvl+1} = [rpos{rlvl+1};trpos{ind(j)}];
        rflux{rlvl+1} = [rflux{rlvl+1};trflux{ind(j)}];
    end
    dlmwrite([ddir,'_rightflux',cname,'.dat'],[rpos{rlvl+1},rflux{rlvl+1}],'precision','%18.14e');
    % Delete data structures
    clear data geometry DoF FE sol flux;
end
% Sum outputs
dlmwrite([ddir,'_fluxsums','.dat'],[dofnum,lfluxsum,rfluxsum,lfluxbndsum,rfluxbndsum],'precision','%18.14e');