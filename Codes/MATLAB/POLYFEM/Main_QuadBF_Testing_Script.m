%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          BF Testing Script
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
clc; close all; format long e; clear persistent; fclose('all');
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Populate global space
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Office');
% Function space
% ------------------------------------------------------------------------------
aa = 1;
bb = 1;
cc = 1;
dd = 1;
ee = 1;
ff = 1;
% Test unit square
% ------------------------------------------------------------------------------
v = [0,0;1,0;1,1;0,1];
vv = [v;.5,0;1,.5;.5,1;0,.5];
f={[1,2],[2,3],[3,4],[4,1]};
n = 101;
[x,y] = meshgrid(linspace(0.000001,.999999,n));
xx=x(:); yy=y(:);
qx = [xx,yy];
[b,g] = mean_value_O2_basis_functions(v,qx,f,2,4);
vvc  = aa*ones(8,1);              qqc  = aa*ones(n*n,1);
vvx  = bb*vv(:,1);                qqx  = bb*qx(:,1);
vvy  = cc*vv(:,2);                qqy  = cc*qx(:,2);
vvx2 = dd*vv(:,1).*vv(:,1);       qqx2 = dd*qx(:,1).*qx(:,1);
vvy2 = ee*vv(:,2).*vv(:,2);       qqy2 = ee*qx(:,2).*qx(:,2);
vvxy = ff*vv(:,1).*vv(:,2);       qqxy = ff*qx(:,1).*qx(:,2);
% Test constant
norm(b*vvc - qqc)
% Test x
norm(b*vvx - qqx)
% Test y
norm(b*vvy - qqy)
% Test x^2
norm(b*vvx2 - qqx2)
% Test y^2
norm(b*vvy2 - qqy2)
% Test xy
norm(b*vvxy - qqxy)
% Test degenerate pentagon
% ------------------------------------------------------------------------------
v = [0,0;1,0;1,1;.5,1;0,1];
vv = [v;.5,0;1,.5;.75,1;.25,1;0,.5];
f={[1,2],[2,3],[3,4],[4,5],[5,1]};
n = 201;
[x,y] = meshgrid(linspace(0,1,n));
xx=x(:); yy=y(:);
qx = [xx,yy];
b = mean_value_O2_basis_functions(v,qx,f,2,5);
vvc  = aa*ones(size(vv,1),1);     qqc  = aa*ones(n*n,1);
vvx  = bb*vv(:,1);                qqx  = bb*qx(:,1);
vvy  = cc*vv(:,2);                qqy  = cc*qx(:,2);
vvx2 = dd*vv(:,1).*vv(:,1);       qqx2 = dd*qx(:,1).*qx(:,1);
vvy2 = ee*vv(:,2).*vv(:,2);       qqy2 = ee*qx(:,2).*qx(:,2);
vvxy = ff*vv(:,1).*vv(:,2);       qqxy = ff*qx(:,1).*qx(:,2);
% Test constant
norm(b*vvc - qqc)
% Test x
norm(b*vvx - qqx)
% Test y
norm(b*vvy - qqy)
% Test x^2
norm(b*vvx2 - qqx2)
% Test y^2
norm(b*vvy2 - qqy2)
% Test xy
norm(b*vvxy - qqxy)