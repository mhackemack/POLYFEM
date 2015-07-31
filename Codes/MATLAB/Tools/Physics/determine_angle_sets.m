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
function data = determine_angle_sets(data, mesh)
% Process some input information
% ------------------------------
dim = mesh.Dimension;
angs = data.Neutronics.Transport.AngularDirections;
% Full Angle Collapsing into a single angle set
% ---------------------------------------------
if strcmp(data.Neutronics.Transport.AngleAggregation, 'all')
    na = data.Neutronics.Transport.NumberAngularDirections;
    data.Neutronics.Transport.AngleSets = {1:na};
    data.Neutronics.Transport.AverageAngles{1} = mean(angs);
    return
end
% No angle collapsing - solve 1 angle at a time (this is the old way)
% -------------------------------------------------------------------
if strcmp(data.Neutronics.Transport.AngleAggregation, 'single')
    na = data.Neutronics.Transport.NumberAngularDirections;
    data.Neutronics.Transport.AngleSets = cell(na, 1);
    for m=1:na
        data.Neutronics.Transport.AngleSets{m} = m;
        data.Neutronics.Transport.AverageAngles(m,:) = angs(m,:);
    end
    return
end
% Switch based on problem dimension
% ---------------------------------
if dim == 1
    ang_sets = determine_1D_angle_sets(data);
elseif dim == 2
    ang_sets = determine_2D_angle_sets(data, mesh);
elseif dim == 3
    ang_sets = determine_3D_angle_sets(data, mesh);
end
% Set Angle Sets
data.Neutronics.Transport.AngleSets = ang_sets;
% Calculate average angle in Angle Set
nangsets = length(ang_sets);
data.Neutronics.Transport.NumberAngleSets = nangsets;
data.Neutronics.Transport.AverageAngles = zeros(nangsets, dim);
for m=1:nangsets
    tangs = ang_sets{m};
    mang = mean(angs(tangs,:),1);
    data.Neutronics.Transport.AverageAngles(m,:) = mang;
end
% Reorder if reflecting boundaries are present
if data.Neutronics.Transport.HasReflectingBoundary
    data = reorder_anglesets(data, mesh);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang_sets = determine_1D_angle_sets(data)
ang_sets = cell(2,1);
% Collect Negative Directions
for m=1:data.Neutronics.Transport.NumberAngularDirections
    if data.Neutronics.Transport.AngularDirections(m,:) < 0
        ang_sets{1} = [ang_sets{1}, m];
    end
end
% Collect Positive Directions
for m=1:data.Neutronics.Transport.NumberAngularDirections
    if data.Neutronics.Transport.AngularDirections(m,:) > 0
        ang_sets{2} = [ang_sets{2}, m];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang_sets = determine_2D_angle_sets(data, mesh)
% NOT Performing Sweep Operations (this is default setting)
if ~data.Neutronics.Transport.performSweeps
    % Collect into 1 angle set if S2
    if data.Neutronics.Transport.NumberAngularDirections == 4
        ang_sets = {1:4};
        return
    end
    % Otherwise, build angle set by direction groupings
    ang_sets = cell(4,1);
    for m=1:data.Neutronics.Transport.NumberAngularDirections
        dir = data.Neutronics.Transport.AngularDirections(m,:);
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
if data.Neutronics.Transport.performSweeps
    % Orthogonal Mesh - Collect into angle chunks
    if mesh.IsOrthogonal
        ang_sets = cell(4,1);
        for m=1:data.Neutronics.Transport.NumberAngularDirections
            dir = data.Neutronics.Transport.AngularDirections(m,:);
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
        ang_sets = cell(data.Neutronics.Transport.NumberAngularDirections, 1);
        for m=1:data.Neutronics.Transport.NumberAngularDirections
            ang_sets{m} = m;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang_sets = determine_3D_angle_sets(data, mesh)
% NOT Performing Sweep Operations (this is default setting)
if ~data.Neutronics.Transport.performSweeps
    % Collect into 1 angle set if S2
    if data.Neutronics.Transport.NumberAngularDirections == 8
        ang_sets = {1:8};
        return
    end
    % Otherwise, build angle set by direction groupings
    ang_sets = cell(8,1);
    for m=1:data.Neutronics.Transport.NumberAngularDirections
        dir = data.Neutronics.Transport.AngularDirections(m,:);
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
if data.Neutronics.Transport.performSweeps
    % Orthogonal Mesh - Collect into angle chunks
    if mesh.IsOrthogonal
        ang_sets = cell(8,1);
        for m=1:data.Neutronics.Transport.NumberAngularDirections
            dir = data.Neutronics.Transport.AngularDirections(m,:);
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
            ang_sets = cell(data.Neutronics.Transport.NumberAngularDirections, 1);
            for m=1:data.Neutronics.Transport.NumberAngularDirections
                ang_sets{m} = m;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = reorder_anglesets(data, mesh)
if data.Neutronics.Transport.HasOpposingReflectingBoundary, return; end
% Reorder anglesets if reflecting boundaries are present
dim = mesh.Dimension;
g_bounds = get_geometry_bounds(mesh);
ref_bounds = data.Neutronics.Transport.ReflectingBoundaries;
ave_dirs = data.Neutronics.Transport.AverageAngles;
nangsets = data.Neutronics.Transport.NumberAngleSets;
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
data.Neutronics.Transport.AngleSets = data.Neutronics.Transport.AngleSets(ind);
data.Neutronics.Transport.AverageAngles = data.Neutronics.Transport.AverageAngles(ind,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_geometry_bounds(mesh)
dim = mesh.Dimension;
out = [mesh.minX, mesh.maxX];
if dim > 1, out = [out;[mesh.minY, mesh.maxY]]; end
if dim > 2, out = [out;[mesh.minZ, mesh.maxZ]]; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%