%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Determine Angle Sets
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = determine_angle_sets_Rev1(data, mesh)
for q=1:length(data.Quadrature)
    data = process_individual_angleset(data, q, mesh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = process_individual_angleset(data, qid, mesh)
% Process some input information
% ------------------------------
dim = mesh.Dimension;
angs = data.Quadrature(qid).AngularDirections;
% Full Angle Collapsing into a single angle set
% ---------------------------------------------
if strcmp(data.Quadrature(qid).AngleAggregation, 'all')
    na = data.Quadrature(qid).NumberAngularDirections;
    data.Quadrature(qid).AngleSets = {1:na};
    data.Quadrature(qid).AverageAngles{1} = mean(angs);
    return
end
% No angle collapsing - solve 1 angle at a time (this is the old way)
% -------------------------------------------------------------------
if strcmp(data.Quadrature(qid).AngleAggregation, 'single')
    na = data.Quadrature(qid).NumberAngularDirections;
    data.Quadrature(qid).AngleSets = cell(na, 1);
    for m=1:na
        data.Quadrature(qid).AngleSets{m} = m;
        data.Quadrature(qid).AverageAngles(m,:) = angs(m,:);
    end
    return
end
% Switch based on problem dimension
% ---------------------------------
if dim == 1
    ang_sets = determine_1D_angle_sets(data, qid);
elseif dim == 2
    ang_sets = determine_2D_angle_sets(data, qid, mesh);
elseif dim == 3
    ang_sets = determine_3D_angle_sets(data, qid, mesh);
end
% Set Angle Sets
data.Quadrature(qid).AngleSets = ang_sets;
% Calculate average angle in Angle Set
nangsets = length(ang_sets);
data.Quadrature(qid).NumberAngleSets = nangsets;
data.Quadrature(qid).AverageAngles = zeros(nangsets, dim);
for m=1:nangsets
    tangs = ang_sets{m};
    mang = mean(angs(tangs,:),1);
    data.Quadrature(qid).AverageAngles(m,:) = mang;
end
% Reorder if reflecting boundaries are present
if data.Transport.HasReflectingBoundary
    data = reorder_anglesets(data, qid, mesh);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang_sets = determine_1D_angle_sets(data, qid)
ang_sets = cell(2,1);
% Collect Negative Directions
for m=1:data.Quadrature(qid).NumberAngularDirections
    if data.Quadrature(qid).AngularDirections(m,:) < 0
        ang_sets{1} = [ang_sets{1}, m];
    end
end
% Collect Positive Directions
for m=1:data.Quadrature(qid).NumberAngularDirections
    if data.Quadrature(qid).AngularDirections(m,:) > 0
        ang_sets{2} = [ang_sets{2}, m];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang_sets = determine_2D_angle_sets(data, qid, mesh)
% NOT Performing Sweep Operations (this is default setting)
if ~data.Transport.PerformSweeps
    % Collect into 1 angle set if S2
    if data.Quadrature(qid).NumberAngularDirections == 4
        ang_sets = {1:4};
        return
    end
    % Otherwise, build angle set by direction groupings
    ang_sets = cell(4,1);
    for m=1:data.Quadrature(qid).NumberAngularDirections
        dir = data.Quadrature(qid).AngularDirections(m,:);
        % (+,+) Directions
        if dir(1) > 0 && dir(2) > 0
            ang_sets{1} = [ang_sets{1}, m];
        end
        % (-,-) Directions
        if dir(1) < 0 && dir(2) < 0
            ang_sets{2} = [ang_sets{2}, m];
        end
        % (+,-) Directions
        if dir(1) > 0 && dir(2) < 0
            ang_sets{3} = [ang_sets{3}, m];
        end
        % (-,+) Directions
        if dir(1) < 0 && dir(2) > 0
            ang_sets{4} = [ang_sets{4}, m];
        end
    end
end
% Performing Sweep Operations
if data.Transport.PerformSweeps
    % Orthogonal Mesh - Collect into angle chunks
    if mesh.IsOrthogonal
        ang_sets = cell(4,1);
        for m=1:data.Quadrature(qid).NumberAngularDirections
            dir = data.Quadrature(qid).AngularDirections(m,:);
            % (+,+) Directions
            if dir(1) > 0 && dir(2) > 0
                ang_sets{1} = [ang_sets{1}, m];
            end
            % (-,-) Directions
            if dir(1) < 0 && dir(2) < 0
                ang_sets{2} = [ang_sets{2}, m];
            end
            % (+,-) Directions
            if dir(1) > 0 && dir(2) < 0
                ang_sets{3} = [ang_sets{3}, m];
            end
            % (-,+) Directions
            if dir(1) < 0 && dir(2) > 0
                ang_sets{4} = [ang_sets{4}, m];
            end
        end
    else
        ang_sets = cell(data.Quadrature(qid).NumberAngularDirections, 1);
        for m=1:data.Quadrature(qid).NumberAngularDirections
            ang_sets{m} = m;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang_sets = determine_3D_angle_sets(data, qid, mesh)
% NOT Performing Sweep Operations (this is default setting)
if ~data.Transport.PerformSweeps
    % Collect into 1 angle set if S2
    if data.Quadrature(qid).NumberAngularDirections == 8
        ang_sets = {1:8};
        return
    end
    % Otherwise, build angle set by direction groupings
    ang_sets = cell(8,1);
    for m=1:data.Quadrature(qid).NumberAngularDirections
        dir = data.Quadrature(qid).AngularDirections(m,:);
        % (+,+,+) Directions
        if dir(1) > 0 && dir(2) > 0 && dir(3) > 0
            ang_sets{1} = [ang_sets{1}, m];
        end
        % (-,-,+) Directions
        if dir(1) < 0 && dir(2) < 0 && dir(3) > 0
            ang_sets{2} = [ang_sets{2}, m];
        end
        % (+,-,+) Directions
        if dir(1) > 0 && dir(2) < 0 && dir(3) > 0
            ang_sets{3} = [ang_sets{3}, m];
        end
        % (-,+,+) Directions
        if dir(1) < 0 && dir(2) > 0 && dir(3) > 0
            ang_sets{4} = [ang_sets{4}, m];
        end
        % (+,+,-) Directions
        if dir(1) > 0 && dir(2) > 0 && dir(3) < 0
            ang_sets{5} = [ang_sets{5}, m];
        end
        % (-,-,-) Directions
        if dir(1) < 0 && dir(2) < 0 && dir(3) < 0
            ang_sets{6} = [ang_sets{6}, m];
        end
        % (+,-,-) Directions
        if dir(1) > 0 && dir(2) < 0 && dir(3) < 0
            ang_sets{7} = [ang_sets{7}, m];
        end
        % (-,+,-) Directions
        if dir(1) < 0 && dir(2) > 0 && dir(3) < 0
            ang_sets{8} = [ang_sets{8}, m];
        end
    end
end
% Performing Sweep Operations
if data.Transport.PerformSweeps
    % Orthogonal Mesh - Collect into angle chunks
    if mesh.IsOrthogonal
        ang_sets = cell(8,1);
        for m=1:data.Quadrature(qid).NumberAngularDirections
            dir = data.Quadrature(qid).AngularDirections(m,:);
            % (+,+,+) Directions
            if dir(1) > 0 && dir(2) > 0 && dir(3) > 0
                ang_sets{1} = [ang_sets{1}, m];
            end
            % (-,-,+) Directions
            if dir(1) < 0 && dir(2) < 0 && dir(3) > 0
                ang_sets{2} = [ang_sets{2}, m];
            end
            % (+,-,+) Directions
            if dir(1) > 0 && dir(2) < 0 && dir(3) > 0
                ang_sets{3} = [ang_sets{3}, m];
            end
            % (-,+,+) Directions
            if dir(1) < 0 && dir(2) > 0 && dir(3) > 0
                ang_sets{4} = [ang_sets{4}, m];
            end
            % (+,+,-) Directions
            if dir(1) > 0 && dir(2) > 0 && dir(3) < 0
                ang_sets{5} = [ang_sets{5}, m];
            end
            % (-,-,-) Directions
            if dir(1) < 0 && dir(2) < 0 && dir(3) < 0
                ang_sets{6} = [ang_sets{6}, m];
            end
            % (+,-,-) Directions
            if dir(1) > 0 && dir(2) < 0 && dir(3) < 0
                ang_sets{7} = [ang_sets{7}, m];
            end
            % (-,+,-) Directions
            if dir(1) < 0 && dir(2) > 0 && dir(3) < 0
                ang_sets{8} = [ang_sets{8}, m];
            end
        end
    else
        if mesh.IsExtruded
            
        else
            ang_sets = cell(data.Quadrature(qid).NumberAngularDirections, 1);
            for m=1:data.Quadrature(qid).NumberAngularDirections
                ang_sets{m} = m;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = reorder_anglesets(data, qid, mesh)
if data.Transport.HasOpposingReflectingBoundary, return; end
% Reorder anglesets if reflecting boundaries are present
dim = mesh.Dimension;
% g_bounds = get_geometry_bounds(mesh);
ref_bounds = data.Transport.ReflectingBoundaries;
ave_dirs = data.Quadrature(qid).AverageAngles;
nangsets = data.Quadrature(qid).NumberAngleSets;
as_order = 1:nangsets;
bound_dirs = zeros(dim, dim, 2);
% Set outward normals
for d=1:dim
    bound_dirs(d,d,1) = -1;
    bound_dirs(d,d,2) = 1;
end
% Loop through angle sets
for m=1:nangsets
    % Loop through dimensions
    for d=1:dim
        for i=1:2
            if ~ref_bounds(d,i), continue; end
            if dot(bound_dirs(d,:,i), ave_dirs(m,:)) < 0
                as_order(m) = as_order(m) + nangsets;
            end
        end
    end
end
% Reorder Anglesets
[~,ind] = sort(as_order);
data.Quadrature(qid).AngleSets = data.Quadrature(qid).AngleSets(ind);
data.Quadrature(qid).AverageAngles = data.Quadrature(qid).AverageAngles(ind,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = get_geometry_bounds(mesh)
% dim = mesh.Dimension;
% out = [mesh.minX, mesh.maxX];
% if dim > 1, out = [out;[mesh.minY, mesh.maxY]]; end
% if dim > 2, out = [out;[mesh.minZ, mesh.maxZ]]; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%