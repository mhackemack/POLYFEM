%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Retrieve the Phase Transformation Matrix
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
function inputs = build_phase_transformation_matrix( data, inputs )
% Preliminary Information
% -----------------------
global glob
dim = data.problem.Dimension;
dim_phase = data.NumberPhasePerDim;
tot_phase = dim_phase^dim;
% Quick Error Checking
% --------------------
% quick_error_checking(dim, dim_phase, inputs.dofs);
% Build Phase Space
wn_norm = 2*pi;
if dim == 1
    p = linspace(0,wn_norm,dim_phase)';
elseif dim == 2
    [px,py] = meshgrid(linspace(0,wn_norm,dim_phase),linspace(0,wn_norm,dim_phase));
    p = [px(:),py(:)];
else
    [px,py,pz] = meshgrid(linspace(0,wn_norm,dim_phase),linspace(0,wn_norm,dim_phase),linspace(0,wn_norm,dim_phase));
    p = [px(:),py(:),pz(:)];
end
% Build all matrices
% ------------------
inputs.phase = cell(inputs.TotalMeshes, 1);
disp('-> Building Phase Matrices.'); rev_str = [];
for m=1:inputs.TotalMeshes
    msg = sprintf('      -> Building Phase for Mesh: %d of %d',m,inputs.TotalMeshes);
    fprintf([rev_str,msg]);
    rev_str = repmat(sprintf('\b'), 1, length(msg));
    
    pt = p;
    gd = get_problem_dimensions(inputs.meshes{m});
    for d=1:dim
        pt(:,d) = pt(:,d)./gd(d,2);
    end
    % Retrieve some mesh/dof information to speed up processing
    % ---------------------------------------------------------
    m_mesh = inputs.meshes{m};
    dof = inputs.dofs{m};
    n_faces = m_mesh.TotalFaces;
    NLocs = dof.NodeLocations;
    FaceID = m_mesh.FaceID;
    FCenter = m_mesh.FaceCenter;
    % Allocate Memory
    inputs.phase{m}.TotalPhases = tot_phase;
    inputs.phase{m}.NumberPhasePerDim = dim_phase;
    inputs.phase{m}.WNList = pt;
    inputs.phase{m}.WN = cell(dim,1);
    inputs.phase{m}.offset = zeros(n_faces, 2);
%     PM  = cell(tot_phase,1);         % this is done to speed up big runs
    % Reform Wave Number Space
    if dim ~= 1
        dd = dim_phase*ones(1,dim);
        for d=1:dim
            inputs.phase{m}.WN{d} = reshape(pt(:,d),dd);
        end
    end
    % Build Phase Matrix Space
    for j=1:tot_phase
%         ptt = pt(j,:)';
%         PM{j} = diag(exp((NLocs*ptt)*1i));
    end
    % Loop through faces
    for f=1:n_faces
        fid = FaceID(f);
        if fid == 0, continue; end
        fcenter = FCenter(f,:);
        % Find periodic dimension
        for d=1:dim
            if sqrt((fcenter(d)-gd(d,1))*(fcenter(d)-gd(d,1))) < glob.small
                pdim = d; gval = -gd(d,2);
                break;
            end
            if sqrt((fcenter(d)-gd(d,2))*(fcenter(d)-gd(d,2))) < glob.small
                pdim = d; gval = gd(d,2);
                break;
            end
        end
        inputs.phase{m}.offset(f,:) = [pdim, gval];
    end
    % Assign Phase Matrices
%     inputs.phase{m}.VolumePM  = PM;
end
fprintf(rev_str); rev_str = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_problem_dimensions( g )
if g.Dimension == 1
    out = [g.minX, g.maxX];
elseif g.Dimension == 2
    out = [g.minX, g.maxX; g.minY, g.maxY];
else
    out = [g.minX, g.maxX; g.minY, g.maxY; g.minZ, g.maxZ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quick_error_checking(dim, n_phase_per_dim, dofs)
global glob
tot_dofs = 0;
for i=1:length(dofs)
    tot_dofs = tot_dofs + dofs{i}.TotalDoFs^2;
end
tot_phase = n_phase_per_dim^dim;
C = 8*3;    % 8 bytes * apparent extra mem storage for structs...
tot_est = C*tot_dofs*tot_phase;
if strcmp(glob.SystemType, 'Home')
    if tot_est > 6e9, error('Total Estimated Counts too high.'); end
elseif strcmp(glob.SystemType, 'Office')
    if tot_est > 2.5e9, error('Total Estimated Counts too high.'); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%