%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Build Fuel Assembly Geometry
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate a geometry object for a nuclear 
%                   reactor assembly.
%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:          All fuel pins will be of the same material. This can be
%                   changed later if desired.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = build_fuel_assembly_geometry(L, num, R, varargin)
% Quick Input Error Checking
num = check_input_data(L, num, R);
data = retrieve_input(varargin{:});
% Build Mesh Components
[nodes, edges, faces] = build_geometry_components(L, num, R, data);
% Pass through Mesh2D Routine
[p, t] = meshfaces(nodes,edges,faces);
% Cleanup Mesh Components
[p, t] = remove_duplicate_geometry(p, t);
% Build Triangulation
tri = triangulation(t,p);
% Build Geometry Object
geometry = GeneralGeometry(2,'Delaunay',tri);
geometry = set_pin_material(geometry, L, num, R, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function num = check_input_data(L, num, R)
num = round(num);
if num < 1, error('Need at least 1 fuel division.'); end
h = L/num;
if 2*R >= h, error('Pin diameter greater than pitch.'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = retrieve_input(varargin)
% Circle Facets
if nargin > 0
    data.num_circ_facets = varargin{1};
else
    data.num_circ_facets = 12;
end
% Maximum Area
if nargin > 1
    data.max_area = varargin{2};
else
    data.max_area = 0.1;
end
% Tolerance
if nargin > 2
    data.tol = varargin{3};
else
    data.tol = 1e-10;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes, edges, faces] = build_geometry_components(L, num, R, data)
% Build Ref Objects
[pin_nodes,pin_edges] = build_ref_pin(R, data.num_circ_facets);
num_pin_tot = size(pin_nodes,1);
% Build Assembly
h = L/num;
faces = cell(num*num+1,1);
nodes = [];
edges = [];
counter = 0; num_pins = 0;
for i=1:num
    x0 = (i-1)/num*L + h/2;
    for j=1:num
        y0 = (j-1)/num*L + h/2;
        % Increment Counter
        counter = counter + 1;
        tpn = [pin_nodes(:,1)+x0, pin_nodes(:,2)+y0];
        tpe = pin_edges + num_pins;
        nodes = [nodes;tpn];
        edges = [edges;tpe];
        faces{counter} = num_pins + (1:num_pin_tot);
        % Increment Pin Numbers
        num_pins = num_pins + num_pin_tot;
    end
end
% Build Bounding Box
tpn = [0,0;L,0;L,L;0,L]; tpe = [1,2;2,3;3,4;4,1];
nodes = [nodes;tpn]; edges = [edges;tpe+size(edges,1)];
faces{end} = 1:size(edges,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges] = build_ref_pin(R, num)
num = num + mod(num,2);
dtheta = pi/(num/2); theta = (-pi:dtheta:(pi-dtheta))';
x = R*cos(theta);
y = R*sin(theta);
nodes = [x, y]; ntot = size(nodes,1);
edges = [(1:ntot)',[(2:ntot)'; 1]];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes, tri] = remove_duplicate_geometry(nodes, tri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = set_pin_material(geometry, L, num, R, data)
% Build Ref Objects
[pin_nodes,~] = build_ref_pin(R, data.num_circ_facets);
num_pin_tot = size(pin_nodes,1);
h = L/num;
for i=1:num
    x0 = (i-1)/num*L + h/2;
    for j=1:num
        y0 = (j-1)/num*L + h/2;
        tpn = [pin_nodes(:,1)+x0, pin_nodes(:,2)+y0];
        geometry.set_cell_matIDs_inside_domain(2,tpn);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%