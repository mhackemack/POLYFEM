%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Split Poly Mesh
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
if exist('pbool', 'var')
    clearvars -except pbool
else
    clear; pbool = false;
end
clc; close all; format long e; clear persistent;
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Begin user input section
% ------------------------------------------------------------------------------
slope = -1;
yintercept = 1;
n = [160000];
% n = [4,16,64,256,1024,4096,16384,65536];
oname = 'SplitPolyMesh_slope=-1_yint=1';
% End user input section
% ------------------------------------------------------------------------------
for g=1:length(n)
    msg = sprintf('Geometry: %d of %d.',g,length(n));
    disp(msg);
    outname = sprintf('geometry_inputs/precompiled/%s_n%d.mat',oname,n(g));
    gname = sprintf('geometry_inputs/precompiled/PolyMesh_SqDomain_L1_n%d.mat',n(g));
    load(gname);
    geometry.split_2d_mesh_on_line(slope,yintercept);
    save(outname,'geometry');
    clear geometry;
end