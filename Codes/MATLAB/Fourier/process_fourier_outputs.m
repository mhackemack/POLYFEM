%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process Fourier Outputs
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
function process_fourier_outputs(data, inputs, outputs)
% Collect some input information
dim      = data.problem.Dimension;
gt       = data.geometry_type;
q_type   = data.Neutronics.Transport.QuadType;
n_sn     = length(data.Neutronics.Transport.SnLevels);
sdm      = data.Neutronics.SpatialMethod;
fdeg     = data.Neutronics.FEMDegree;
t_type   = [data.Neutronics.TransportMethod,'_',data.Neutronics.Transport.transportType];
dsa_type = data.Neutronics.DSAType;
if data.Neutronics.FEMLumping
    lump = 'L';
else
    lump = 'U';
end
% Adjust IP/MIP Name with IP coefficient
if strcmp(dsa_type,'MIP') ||strcmp(dsa_type,'IP')
    dsa_type = [dsa_type,'_C=',num2str(data.Neutronics.IP_Constant)];
end
f_name = [t_type,'_',dsa_type,'_'];
sigt = data.Neutronics.Transport.TotalXS;
mfp = inputs.x * sigt;
% Preliminary plotting setup
if data.Output.plotting_bool
    figure(1); hFig = figure(1); hold on;
    set(hFig,'Position',[1,1,1200,700]);
    plot_counter = 0;
end
% Loop through different quadrature levels
for m=1:n_sn
    counter = 0;
    mm = data.Neutronics.Transport.SnLevels(m);
    quad_name = [q_type, num2str(mm)];
    % Dimension = 1
    if dim == 1
        dn = 'outputs/1D/';
        te = zeros(inputs.nx,1);
        for i=1:inputs.nx
            counter = counter + 1;
            te(i) = outputs{m,counter}.Eigen.Max;
        end
        % Plot here
        if data.Output.plotting_bool
            plot_counter = plot_counter + 1;
            plot(mfp, te);
            leg_names{plot_counter} = ['S',num2str(mm)];
        end
        % Print here
        fname = [dn,f_name,'S',num2str(mm)];
        if data.Output.file_bool, dlmwrite([fname,'.dat'],[mfp,te],'precision','%14.8e'); end
    % Dimension = 2
    elseif dim == 2
        dn = ['outputs/2D/',gt,'/'];
        % Loop through X/Y aspect ratios
        for j=1:inputs.nyz
            te = zeros(inputs.nx,1);
            % Loop through meshes
            for i=1:inputs.nx
                counter = counter + 1;
                te(i) = outputs{m,counter}.Eigen.Max;
            end
            % Plot here
            if data.Output.plotting_bool
                plot_counter = plot_counter + 1;
                plot(mfp, te);
                leg_names{plot_counter} = ['S',num2str(mm),', X/Y = ',num2str(inputs.yz(j))];
            end
            % Print here
            fname = [dn,f_name,quad_name,'_',lump,sdm,'_k',num2str(fdeg),'_XY=',num2str(inputs.yz(j))];
            if data.Output.file_bool, dlmwrite([fname,'.dat'],[mfp,te],'precision','%14.8e'); end
        end
    % Dimension = 3
    elseif dim == 3
        dn = ['outputs/3D/',gt,'/'];
        % Loop through X/YZ aspect ratios
        for j=1:inputs.nyz
            te = zeros(inputs.nx,1);
            % Loop through meshes
            for i=1:inputs.nx
                counter = counter + 1;
                te(i) = outputs{m,counter}.Eigen.Max;
            end
            % Plot here
            if data.Output.plotting_bool
                plot_counter = plot_counter + 1;
                plot(mfp, te);
                leg_names{plot_counter} = ['S',num2str(mm),', X/YZ = ',num2str(inputs.yz(j))];
            end
            % Print here
            fname = [dn,f_name,quad_name,'_',lump,sdm,'_k',num2str(fdeg),'_XYZ=',num2str(inputs.yz(j))];
            if data.Output.file_bool, dlmwrite([fname,'.dat'],[mfp,te],'precision','%14.8e'); end
        end
    end
end


% Fixup plot properties
% ---------------------
if data.Output.plotting_bool
    box on; font_size = 22;
    set(hFig,'Position',[1,1,1200,700])
    set(gca,'xscale','log');
    set(gca,'XGrid','on','XMinorGrid','off');
    set(gca,'YGrid','on','YMinorGrid','off');
    set(gca,'FontName','Times New Roman','FontSize',font_size);
    xlabel('Cell Size in MFP', 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
    ylabel('Spectral Radius', 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
end





