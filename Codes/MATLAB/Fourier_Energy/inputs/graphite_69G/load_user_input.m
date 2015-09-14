%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Graphite 69 Energy Groups
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
function data = load_user_input()
% Get Cross Section Data
% ------------------------------------------------------------------------------
% Total XS
load('MT_1.mat');
data.TotalXS = mat;
% Scattering XS
load('MT_2500.mat');
data.ScatteringXS = mat;
% Energy Bounds
load('EnergyBounds.mat');
data.EnergyBounds = mat;
% Define Quadrature Set
% ------------------------------------------------------------------------------
data.Neutronics.Transport.fluxMoments = 0;
data.Neutronics.Transport.AngleAggregation = 'auto';
data.Neutronics.Transport.QuadType = 'LS';
data.Neutronics.Transport.SnLevels = 8;
data.Neutronics.Transport.PolarLevels = 4;
data.Neutronics.Transport.AzimuthalLevels = 4;
% Define Energy and Iteration Structure
% ------------------------------------------------------------------------------

