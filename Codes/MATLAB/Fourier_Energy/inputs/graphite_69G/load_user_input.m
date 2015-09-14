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
% Input/Output Information
% ------------------------------------------------------------------------------
data.IO.OutputDirectory = 'graphite_69G';
% Get Cross Section Data
% ------------------------------------------------------------------------------
% Total XS
load('MT_1.mat');
data.XS.TotalXS = mat;
% Scattering XS
load('MT_2500.mat');
data.XS.ScatteringXS = mat;
% Energy Bounds
load('EnergyBounds.mat');
data.XS.EnergyBounds = mat;
% Define Quadrature Set
% ------------------------------------------------------------------------------
data.Quad.fluxMoments = 0;
data.Quad.AngleAggregation = 'auto';
data.Quad.QuadType = 'LS';
data.Quad.SnLevels = 8;
data.Quad.PolarLevels = 4;
data.Quad.AzimuthalLevels = 4;
% Define Energy and Iteration Structure
% ------------------------------------------------------------------------------

