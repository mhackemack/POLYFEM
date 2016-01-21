%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Add XS Component
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
function data = add_xs_component_contribution(data,XSID,matid,xs_name,density)
global glob
% Quick Input Checking
% ------------------------------------------------------------------------------
if nargin < 4, error('Not engough input arguments.'); end
if nargin < 5, density = 1.0; end
nm   = data.problem.NumberMaterials;
ng   = data.Groups.NumberEnergyGroups;
nmom = data.Transport.PnOrder+1;
% Check if XS field has been built
% ------------------------------------------------------------------------------
if ~isfield(data, 'XS'), data.XS = make_empty_XS_field(nm,ng,nmom); end
if length(data.XS) < XSID
    for i=length(data.XS)+1:XSID, data.XS(i) = make_empty_XS_field(nm,ng,nmom); end
end
% Check if XS Files exist
% ------------------------------------------------------------------------------
xs_dir = [glob.xs_dir, xs_name];
if ~exist(xs_dir, 'dir'), error('XS component does not exist.'); end
% Add XS contribution
% ------------------------------------------------------------------------------
data.XS(XSID) = add_comp_contribution(data.XS(XSID),xs_dir,matid,density);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = make_empty_XS_field(nm,ng,nmom)
out.TotalXS = zeros(nm,ng);
out.AbsorbXS = zeros(nm,ng);
out.ScatteringXS = zeros(nm,ng,ng,nmom);
out.FissionXS = zeros(nm,ng);
out.NuBar = zeros(nm,ng);
out.FissSpec = zeros(nm,ng);
out.ExtSource = [];
out.BCFlags = [];
out.BCVals = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XS = add_comp_contribution(XS, xs_dir, matid, density)
% Total XS contribution - MT1
if exist([xs_dir,'/MT_1.mat'], 'file')
    T=retrieve_xs_component_from_file([xs_dir,'/MT_1.mat']);
    if size(T,1) > size(T,2), T = T'; end
    XS.TotalXS(matid,:) = XS.TotalXS(matid,:) + density*T;
end
% Scattering XS contribution
nsxs = size(XS.ScatteringXS,4);
if exist([xs_dir,'/MT_2500.mat'], 'file')
    S=retrieve_xs_component_from_file([xs_dir,'/MT_2500.mat']);
    XS.ScatteringXS(matid,:,:,:) = squeeze(XS.ScatteringXS(matid,:,:,:)) + density*S(:,:,1:nsxs);
elseif exist([xs_dir,'/MT_2501.mat'], 'file')
    S=retrieve_xs_component_from_file([xs_dir,'/MT_2501.mat']);
    XS.ScatteringXS(matid,:,:,:) = squeeze(XS.ScatteringXS(matid,:,:,:)) + density*S(:,:,1:nsxs);
end
% Absorption XS contribution - check if MT27 exists, otherwise build it from the
% the total and scattering cross sections
if exist([xs_dir,'/MT_27.mat'], 'file')
    T=retrieve_xs_component_from_file([xs_dir,'/MT_27.mat']);
    if size(T,1) > size(T,2), T = T'; end
    XS.AbsorbXS(matid,:) = XS.AbsorbXS(matid,:) + density*T;
else
    
end
% Fission XS - MT18
if exist([xs_dir,'/MT_18.mat'], 'file')
    T=retrieve_xs_component_from_file([xs_dir,'/MT_18.mat']);
    if size(T,1) > size(T,2), T = T'; end
    XS.FissionXS(matid,:) = XS.FissionXS(matid,:) + density*T;
end
% Nu-Bar - MT452
if exist([xs_dir,'/MT_452.mat'], 'file')
    T=retrieve_xs_component_from_file([xs_dir,'/MT_452.mat']);
    if size(T,1) > size(T,2), T = T'; end
    XS.NuBar(matid,:) = XS.NuBar(matid,:) + density*T;
end
% Fission spectrum (Chi) - MT1018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = retrieve_xs_component_from_file(xs_file_name)
M = open(xs_file_name);
if isnumeric(M)
    out = M;
elseif isstruct(M)
    if isfield(M,'mat'), out = M.mat;
    else error('Cannot determine xs field.');
    end
elseif isempty(M)
    error('No data is present in the xs file.');
else
    error('Cannot determine xs file organization.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%