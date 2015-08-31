%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          PDT XS Reader (Script)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    MATLAB script to read PDT XS files and pull pertinent 
%                   data out. The main XS of interest are the total XS and 
%                   the scattering matrices.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear; clc; close all; format long e
% Set all user input information necessary
% ------------------------------------------------------------------------------
data.file_name = 'XSFiles/pdt-119g-graphite_cnat.cx';
data.out_dir = 'inputs/graphite_119G';
data.num_groups = 119;          % removes some parsing burden
data.iscat = 8;                 % make sure this one is correct
data.scatt_enums = 2500;        % scattering kernels
data.enums_1G = [1099,1,2,4];   % 1G XS numbers
% ------------------------------------------------------------------------------
% Begin program execution
% print_XSR_heading();
% f_out = xs_parse_file(data.file_name);
xs_data = xs_strip_data(data, f_out);
% xs_output_data(data, xs_data);
% ------------------------------------------------------------------------------
