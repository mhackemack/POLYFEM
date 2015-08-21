%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          DG Geometry Generation Script
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
clear;clc;close all;
fclose('all');

in_dir = 'geometry_inputs/dg_raw/';
out_dir = 'geometry_inputs/precompiled/';

files = dir( fullfile(in_dir,'*.txt') );
fsizes = {files.bytes}';
files = {files.name}';

for i=1:numel(files)
    [pathstr,name,ext] = fileparts([in_dir,files{i}]);
    g_size = fsizes{i};
    g_time = 1.09776E-10*g_size*g_size + 8.16314E-06*g_size + 2.17143E-03;
    disp(['Generating File #:  ',num2str(i)])
    disp(['File Name:          ', name])
    disp(['File Size in Bytes: ', num2str(g_size)])
    disp(['Estimated Time:     ',num2str(g_time)])
    disp(' ')
    geometry = GeneralGeometry(2,'DG',[in_dir,files{i}]);
    save([out_dir,name,'.mat'], 'geometry');
    clear geometry
end
