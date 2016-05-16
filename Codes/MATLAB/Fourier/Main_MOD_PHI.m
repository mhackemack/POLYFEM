%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Fourier Analysis Script - Modified PHI Problem
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
inp = '2D_MOD_PHI'; addpath([glob.input_path,inp]);
% Problem inputs
% sigt = [10000,100000,1000000];
% c = 0.9999;
sigt = [10];
c    = [0.9,0.99,0.999,0.9999,0.99999,0.999999];
ncells = 2;
L = 1;
ngrid = 101;
data = load_user_input(ncells,L);
% end user input section
% ------------------------------------------------------------------------------
% Compute all solutions
% ---------------------
print_FA_heading(data);
b_func = get_build_function(data);
outdir = sprintf('outputs/MODPHI/%dD/%s/',data.problem.Dimension,data.geometry.type);
if data.Neutronics.FEMLumping
    lump = 'L';
else
    lump = 'U';
end
if ~data.Neutronics.PerformAcceleration % Unaccelerated
    outname = sprintf('%s_%s%s%d_ncells=%d_L=%d',data.Neutronics.TransportMethod,lump,data.Neutronics.SpatialMethod,data.Neutronics.FEMDegree,ncells,L);
elseif data.Neutronics.PerformAcceleration % Accelerated
    outname = sprintf('%s_%s_C=%d_%s%s%d_ncells=%d_L=%d',data.Neutronics.TransportMethod,data.Neutronics.DSAType,data.Neutronics.IP_Constant,lump,data.Neutronics.SpatialMethod,data.Neutronics.FEMDegree,ncells,L);
end
if ~isequal(exist(outdir, 'dir'),7),mkdir(outdir); end
SR = zeros(length(sigt),length(c),length(data.Neutronics.Transport.SnLevels));
% Loop through total cross sections
for t=1:length(sigt)
    data = load_user_input(ncells,L);
    for cc=1:length(c)
        data = set_phi_xs(data, sigt(t), c(cc));
        % Build input space
        [data, inputs] = process_fourier_inputs( data );
        % Run the Search problems
        data.Type = 'Search';
        data.NumberPhasePerDim = 5;
        pmin = sqrt(eps); pmax = 2*pi - sqrt(eps);
        data.PhaseXSpacing = linspace(pmin,pmax,data.NumberPhasePerDim);
        data.PhaseYSpacing = linspace(pmin,pmax,data.NumberPhasePerDim);
        inputs = build_phase_transformation_matrix(data, inputs);
        outputs = calculate_eigenspectrums(data, inputs);
        % Process Search outputs - also calculate maximum case for eigenspectrum
        for q=1:length(data.Neutronics.Transport.SnLevels)
            qlvl = data.Neutronics.Transport.SnLevels(q);
            % Save maximum eigenvalue - compute largest eigenvalue distribution
            SR(t,cc,q) = outputs{q,1}.Eigen.Max;
            [inp, ~] = combine_input_set(data, inputs, 1, q);
            P = b_func(outputs{q,1}.Eigen.MaxLambda,inp);
            [V,D] = eig(P); D=diag(D);
            % Save output data
            fulloutname = sprintf('%s_%s%d_sigt=%d_c=%g',outname,data.Neutronics.Transport.QuadType,qlvl,sigt(t),c(cc));
            dlmwrite([outdir,fulloutname,'_Eigenvalues.dat'],[D,real(D),imag(D)]);
            dlmwrite([outdir,fulloutname,'_Eigenvectors.dat'],V);
        end
%         % Run the Grid problems
%         data.Type = 'Grid';
%         data.NumberPhasePerDim = ngrid;
%         pmin = 0; pmax = 2*pi;
%         data.PhaseXSpacing = linspace(pmin,pmax,data.NumberPhasePerDim);
%         data.PhaseYSpacing = linspace(pmin,pmax,data.NumberPhasePerDim);
%         inputs = build_phase_transformation_matrix(data, inputs);
%         outputs = calculate_eigenspectrums(data, inputs);
%         % Process Grid outputs
%         for q=1:length(data.Neutronics.Transport.SnLevels)
%             qlvl = data.Neutronics.Transport.SnLevels(q);
%             % Make contour plot
%             x = inputs.phase{1}.WN{1}; x = x./max(max(x))*2*pi;
%             y = inputs.phase{1}.WN{2}; y = y./max(max(y))*2*pi;
%             z = outputs{q,1}.Eigen.Grid;
%             contourf(x,y,z,20); colorbar;
%             xlabel('\lambda_x');
%             ylabel('\lambda_y');
%             set(gca,'FontSize',11);
% %             tightfig(gcf); colorbar;
%             % Save contour plot and Grid data
%             fulloutname = sprintf('%s_%s%d_sigt=%d_c=%g',outname,data.Neutronics.Transport.QuadType,qlvl,sigt(t),c(cc));
%             dlmwrite([outdir,fulloutname,'_GridData',num2str(data.NumberPhasePerDim),'.dat'],z);
%             savefig(gcf,[outdir,fulloutname,'_contour.fig']);
%             print(gcf,'-dpng',[outdir,fulloutname,'_contour.png']);
% %             print(gcf,'-depsc',[outdir,fulloutname,'_contour.eps']);
%             close(gcf);
%         end
    end
end
% Save off grid information
% ------------------------------------------------------------------------------
x = inputs.phase{1}.WN{1}; x = x./max(max(x))*2*pi;
y = inputs.phase{1}.WN{2}; y = y./max(max(y))*2*pi;
dlmwrite([outdir,'GridX',num2str(ngrid),'.dat'],x);
dlmwrite([outdir,'GridY',num2str(ngrid),'.dat'],y);
% % Save off maximum spectral radii information
% % ------------------------------------------------------------------------------
for q=1:length(data.Neutronics.Transport.SnLevels)
    qlvl = data.Neutronics.Transport.SnLevels(q);
    fulloutname = sprintf('%s_%s%d',outname,data.Neutronics.Transport.QuadType,qlvl);
    dlmwrite([outdir,fulloutname,'_SpectralRadii.dat'],[0,c;sigt',SR(:,:,q)]);
end

