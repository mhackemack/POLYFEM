%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get Reference Polygons
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    Builds the 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = get_reference_polygons(vlen)
% Get total number of polygons to build
num_polys = length(vlen);
% Get polygon arrays
% ------------------------------------------------------------------------------
verts = cell(num_polys,1); faces = cell(num_polys,1);
% Loop through polygons
for i=1:num_polys
    % Skip if number of verts is less than 3.
    tvlen = vlen(i);
    if tvlen < 3, continue; end
    % Get some geometric properties
    int_angle = (180 - (tvlen-2)*180/tvlen)*(pi/180);
    % Build polygon geometric structures
    if tvlen == 3
        verts{i} = [0,0;1,0;0,1];
        faces{i} = {[1,2],[2,3],[3,1]};
    elseif tvlen == 4
        verts{i} = [0,0;1,0;1,1;0,1];
        faces{i} = {[1,2],[2,3],[3,4],[4,1]};
    elseif tvlen > 4
        verts{i} = zeros(tvlen,2);
        faces{i} = cell(1,tvlen);
        % Loop through vertex numbers
        tang = pi/2; % starts the vertices at the top
        for j=1:tvlen
            faces{i}{j} = [j,mod(j,tvlen)+1];
            verts{i}(j,:) = [cos(tang), sin(tang)];
            tang = tang + int_angle;
        end
    end
end
% Set output arguments
% ------------------------------------------------------------------------------
varargout{1} = verts;
varargout{2} = faces;
