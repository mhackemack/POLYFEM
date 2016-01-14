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
% Get some problem 

% Check if XS field has been built
% ------------------------------------------------------------------------------
if ~isfield(data, 'XS'), data.XS = make_empty_XS_field(); end
if length(data.XS) < XSID
    for i=length(data.XS)+1:XSID, data.XS(i) = make_empty_XS_field(); end
end
% Check if XS Files exist
% ------------------------------------------------------------------------------
xs_dir = [glob.xs_dir, xs_name];

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XS = add_comp_contribution(XS, xs_dir, density)
% Total XS contribution - MT1
if ~exist([xs_dir,'MT_1.mat'], 'file'), error(''); end
% Scattering XS contribution

% Absorption XS contribution - check if MT27 exists, otherwise build it from the
% the total and scattering cross sections
if ~exist([xs_dir,'MT_27.mat'], 'file')
    
else
    
end
% Fission XS - MT18

% Nu-Bar - MT452

% Fission spectrum (Chi) - MT1018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%