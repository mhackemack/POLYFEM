%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          BF Timing Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
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
clc; close all; format long e;
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Populate global space
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Office');
% Being User Input Section
% ------------------------------------------------------------------------------
outdir = 'outputs/BF_Times/';
names  = {'Wachspress','PWL','MV','ME'};
funcs  = {@wachspress_basis_functions,@PWLD_basis_functions,@mean_value_basis_functions,@max_entropy_basis_functions};
fedegs = [1,2];
pnums  = 3:12;
qnums  = 1:8;
num_iters = 100;
% Allocate Memory Space
% ------------------------------------------------------------------------------
numqpts = zeros(length(pnums), length(qnums));
meantimes = zeros(length(pnums), length(qnums), length(fedegs), length(names));
stdtimes = zeros(length(pnums), length(qnums), length(fedegs), length(names));
% Run Problem Suite
% ------------------------------------------------------------------------------
% Loop through polygons
for p=1:length(pnums)
    fprintf('Polygon %d of %d.\n',p,length(pnums))
    [v, f] = RegularPolygon(pnums(p),1);
    % Loop through quadrature orders
    for q=1:length(qnums)
        fprintf('  Quadrature %d of %d.\n',q,length(qnums))
        % Retrieve volume quadrature
        [qx, qw] = get_general_volume_quadrature(v,f,qnums(q),true);
        numqpts(p,q) = length(qw);
        % Loop through fem degrees
        for k=1:length(fedegs)
            % Loop through basis functions
            for b=1:length(names)
                fprintf('    Degree %d of %d, BF %d of %d.\n',k,length(fedegs),b,length(names))
                % Allocate temporary timing arrays
                ttimes = zeros(num_iters,1);
                % Loop through iterations
                for i=1:num_iters
                    t = tic;
                    funcs{b}(v,qx,f,fedegs(k),pnums(p));
                    ttimes(i) = toc(t);
                end
                % Calcualte standard deviation
                meantimes(p,q,k,b) = mean(ttimes);
                stdtimes(p,q,k,b) = sqrt(sum((ttimes - meantimes(p,q,k,b)).^2)/num_iters);
            end
        end
    end
end
% Save off timing data
% ------------------------------------------------------------------------------
dlmwrite(sprintf('%sNumberQuadraturePoints.dat',outdir),[0,qnums;pnums',numqpts]);
% Loop through fem degrees
for k=1:length(fedegs)
    % Loop through basis functions
    for b=1:length(names)
        dlmwrite(sprintf('%s%s_k%d_meantimes.dat',outdir,names{b},fedegs(k)), meantimes(:,:,k,b),'precision','%10.8e');
        dlmwrite(sprintf('%s%s_k%d_std.dat',outdir,names{b},fedegs(k)), stdtimes(:,:,k,b),'precision','%10.8e');
    end
end