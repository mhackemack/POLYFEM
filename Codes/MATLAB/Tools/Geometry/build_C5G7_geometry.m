%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Build C5G7 Benchmark Geometry
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to generate a geometry object for the C5G7
%                   benchmark case.
%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geometry = build_C5G7_geometry(varargin)
% Quick Input Checking
data = retrieve_input(varargin{:});
% Build Mesh Components
[nodes, edges, faces] = build_geometry_components(data);
% Pass through Mesh2D Routine
[p, t] = meshfaces(nodes,edges,faces);
% Cleanup Mesh Components
[p, t] = remove_duplicate_geometry(p, t);
% Build Triangulation
tri = triangulation(t,p);
% Build Geometry Object
geometry = GeneralGeometry(2,'Delaunay',tri);

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
function [nodes, edges, faces] = build_geometry_components(data)
R = 0.54; h = 1.26;
LL = 21.42; L = 64.26;
% Build Ref Objects
[pin_nodes,pin_edges] = build_ref_pin(R, data.num_circ_facets);
[ass_nodes,ass_edges,ass_faces] = build_ref_assembly(pin_nodes, pin_edges);
% Combine all assemblies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges] = build_ref_pin(R, num)
num = num + mod(num,2);
dtheta = pi/(num/2); theta = (-pi:dtheta:(pi-dtheta))';
x = R*cos(theta);
y = R*sin(theta);
nodes = [x, y]; ntot = size(nodes,1);
edges = [(1:ntot)',[(2:ntot)'; 1]];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes, edges, faces] = build_ref_assembly(pin_nodes, pin_edges)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes, tri] = remove_duplicate_geometry(nodes, tri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%