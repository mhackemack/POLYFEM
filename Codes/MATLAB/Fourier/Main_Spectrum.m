%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Fourier Analysis Script - Search Method
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
clc; close all; format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Define Path
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Home');
glob.print_info = false;
% Load all user inputs
% ------------------------------------------------------------------------------
inp = '1D_99G_WGS_DSA'; addpath([glob.input_path,inp]);
data = load_user_input();
% additional inputs
data.Type = 'Search';
data.NumberPhasePerDim = 3;
% data.PhaseXSpacing = linspace(0,.5,121);
% data.PhaseYSpacing = linspace(0,.5,121);
% end user input section
% ------------------------------------------------------------------------------
% Populate data and output structures
% -----------------------------------
print_FA_heading(data);
[data, inputs] = process_fourier_inputs( data );
inputs = build_phase_transformation_matrix(data, inputs);
% Retrieve all spectrum data and postprocess
% ------------------------------------------------------------------------------
b_func = get_build_function(data);
outputs = calculate_eigenspectrums(data, inputs);
Dout = cell(length(data.Neutronics.Transport.SnLevels),inputs.TotalMeshes);
Vout = cell(length(data.Neutronics.Transport.SnLevels),inputs.TotalMeshes);
Dmax = zeros(length(data.Neutronics.Transport.SnLevels),inputs.TotalMeshes);
for q=1:length(data.Neutronics.Transport.SnLevels)
    for m=1:inputs.TotalMeshes
        [inp, ~] = combine_input_set(data, inputs, m, q);
        P = b_func(outputs{q,m}.Eigen.MaxLambda,inp);
        [V,D] = eig(P); D=diag(D);
        Dout{q,m} = D;
        Vout{q,m} = V;
        [Dmax(q,m),ind] = max(abs(D));
        Vmax(:,q,m) = V(:,ind);
%         Dmeshes = [Dmeshes,real(D),imag(D)];
        % Plot data if specified in input
        if data.Output.plotting_bool
            
        end
        % Output data to file if specified in input
        if data.Output.file_bool
            
        end
    end
end
% ------------------------------------------------------------------------------