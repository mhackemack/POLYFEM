%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Reference Polygon Disk Mapper Script
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
clc; close all; format long e; clear;
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Office');
% Being User Input Section
% ------------------------------------------------------------------------------
polynums = 3:30;
sdm = {'PWLD','WACHSPRESS','MV','MAXENT'};
fedeg = [1,2];
quadoffset = 0:6;
% ---
outdir = 'SCCM';
% Execute suite and generate DiskMapper classes
% ------------------------------------------------------------------------------
% Loop through polygons
for p=1:length(polynums)
    pnum = polynums(p);
    % Loop through basis functions
    for b=1:length(sdm)
        bfname = sdm{b};
        % Loop through basis function orders
        for k=1:length(fedeg)
            k2 = 2*fedeg(k);
            % Loop through quadrature orders
            for q=1:length(quadoffset)
                qq = k2 + quadoffset(q);
                % Build Reference Polygon Disk Mapper object
                RPDM = RefPolygonDiskMapper(pnum,bfname,k,qq);
            end
        end
    end
end
            