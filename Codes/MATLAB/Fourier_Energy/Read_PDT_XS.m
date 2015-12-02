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
clear; clc; close all; format long e
% Set all user input information necessary
% ------------------------------------------------------------------------------
data.file_name = 'Z:\PDT\pdt_inputs\work_inputs\Simple_TG_Testing\pdt-99g_ni58.cx';
data.out_dir = 'inputs/Ni58_99G';
data.num_groups = 99;           % removes some parsing burden
data.iscat = 8;                 % make sure this one is correct
data.scatt_enums = [2500,2501,2519]; % scattering kernels
data.enums_1G = [1099,1,2,4];   % 1G XS numbers
% ------------------------------------------------------------------------------
% Begin program execution
print_XSR_heading();
f_out = xs_parse_file(data.file_name);
xs_data = xs_strip_data(data, f_out);
xs_output_data(data, xs_data);
% ------------------------------------------------------------------------------
