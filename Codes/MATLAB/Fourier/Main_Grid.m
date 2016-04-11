%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Fourier Analysis Script
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
clc; close all; %format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Define Path
% -----------
global glob
glob = get_globals('Office');
glob.print_info = false;
% Define all user inputs
% ------------------------------------------------------------------------------
inp = '2D_1G_DSA'; addpath([glob.input_path,inp]);
data = load_user_input();
% additional inputs
data.Type = 'Grid';
n = 2;
data.NumberPhasePerDim = n;
% pmin = 0; pmax = 1/10;
pmin = 0; pmax = 2*pi;
% wn_norm = 2*pi; pmin = sqrt(eps); pmax = wn_norm - sqrt(eps);
data.PhaseXSpacing = linspace(pmin,pmax,data.NumberPhasePerDim);
data.PhaseYSpacing = linspace(pmin,pmax,data.NumberPhasePerDim);
% end user input section
% ------------------------------------------------------------------------------
% Populate data and output structures
% -----------------------------------
print_FA_heading(data);
[data, inputs] = process_fourier_inputs( data );
inputs = build_phase_transformation_matrix(data, inputs);
% Create directory and output file names
% --------------------------------------
outdir = sprintf('outputs/Grid/%dD/%s/',data.problem.Dimension,data.geometry.type);
if data.Neutronics.FEMLumping
    lump = 'L';
else
    lump = 'U';
end
if ~data.Neutronics.PerformAcceleration % Unaccelerated
    outname = sprintf('%s_%s',data.Neutronics.TransportMethod,lump,data.Neutronics.SpatialMethod);
elseif data.Neutronics.PerformAcceleration % Accelerated
    outname = sprintf('%s_%s_C=%d_%s%s%d',data.Neutronics.TransportMethod,data.Neutronics.DSAType,data.Neutronics.IP_Constant,lump,data.Neutronics.SpatialMethod,data.Neutronics.FEMDegree);
end
if ~isequal(exist(outdir, 'dir'),7),mkdir(outdir); end
% Retrieve all spectrum data and postprocess
% ------------------------------------------------------------------------------
outputs = calculate_eigenspectrums(data, inputs);
% Loop through quadrature sets
if data.Output.plotting_bool
    for q=1:length(data.Neutronics.Transport.SnLevels)
        qlvl = data.Neutronics.Transport.SnLevels(q);
        % Loop through meshes in 1D
        if data.problem.Dimension == 1
            for i=1:TotalMeshes
                
            end
        elseif data.problem.Dimension == 2
            % Loop through meshes in 2D
            c = 0;
            for j=1:inputs.nyz
                yy = inputs.yz(j);
                for i=1:inputs.nx
                    xx = inputs.x(i);
                    c = c + 1;
                    x = inputs.phase{c}.WN{1}; x = x./max(max(x))*(max(data.PhaseXSpacing));
                    y = inputs.phase{c}.WN{2}; y = y./max(max(y))*(max(data.PhaseYSpacing));
                    v = outputs{q,c}.Eigen.Grid;
                    contourf(x,y,v,20); colorbar;
                    xlabel('\lambda_x');
                    ylabel('\lambda_y');
                    set(gca,'FontSize',11);
                    % Save contour plot and Grid data
                    fulloutname = sprintf('%s_%s%d_x=%g_dydx=%g',outname,data.Neutronics.Transport.QuadType,qlvl,xx,yy);
                    dlmwrite([outdir,fulloutname,'_GridData',num2str(data.NumberPhasePerDim),'.dat'],v);
                    savefig(gcf,[outdir,fulloutname,'_contour.fig']);
                    print(gcf,'-dpng',[outdir,fulloutname,'_contour.png']);
                    print(gcf,'-depsc',[outdir,fulloutname,'_contour.eps']);
                    close(gcf);
                end
            end
        elseif data.problem.Dimension == 3
            % Loop through meshes in 3D
            c = 0;
            for j=1:inputs.nyz
                for i=1:inputs.nx
                    c = c + 1;
                    x = inputs.phase{c}.WN{1};
                    y = inputs.phase{c}.WN{2};
                    z = inputs.phase{c}.WN{3};
                    v = outputs{q,c}.Eigen.Grid;
                end
            end
        end
    end
end
