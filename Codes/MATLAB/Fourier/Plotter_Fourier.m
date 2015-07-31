% Clear Project Space
% -------------------
if exist('pbool', 'var')
    clearvars -except pbool
else
    clear; pbool = false;
end
clc; close all; format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Populate global space
% ---------------------
global glob
glob = get_globals('Office');
dir_root_name = 'outputs/';
out_name = [dir_root_name,'figures/'];
% Begin User Input Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_sets = 4; p_info = cell(tot_sets, 1);
save_name = 'SI_MIP_hex_C=1_LS2,4,8,16_y=1';
legend_loc = 'northwest'; font_size = 22;
yaxis_lim = 'manual'; yaxis_max = 0.7;
xaxis_lim = 'auto'; xaxis_b = [1e-2,1e4];
% 1st Info Set
p_info{1} = struct();
p_info{1}.dim = 3;
p_info{1}.g_type = 'cart';
p_info{1}.name = 'SI_MIP_C=1_LS2_PWLD_XYZ=1';
p_info{1}.title = 'S2';
p_info{1}.symbol = 'k-';
% 2nd Info Set
p_info{2} = struct();
p_info{2}.dim = 3;
p_info{2}.g_type = 'cart';
p_info{2}.name = 'SI_MIP_C=1_LS4_PWLD_XYZ=1';
p_info{2}.title = 'S4';
p_info{2}.symbol = 'b-';
% 3rd Info Set
p_info{3} = struct();
p_info{3}.dim = 3;
p_info{3}.g_type = 'cart';
p_info{3}.name = 'SI_MIP_C=1_LS8_PWLD_XYZ=1';
p_info{3}.title = 'S8';
p_info{3}.symbol = 'r-';
% 4th Info Set
p_info{4} = struct();
p_info{4}.dim = 3;
p_info{4}.g_type = 'cart';
p_info{4}.name = 'SI_MIP_C=1_LS16_PWLD_XYZ=1';
p_info{4}.title = 'S16';
p_info{4}.symbol = 'm-';
% 5th Info Set
% p_info{5} = struct();
% p_info{5}.dim = 3;
% p_info{5}.g_type = 'cart';
% p_info{5}.name = 'SI_MIP_C=4_LS8_PWLD_XYZ=4';
% p_info{5}.title = 'Y/X = 4';
% p_info{5}.symbol = 'r--';
% 6th Info Set
% p_info{6} = struct();
% p_info{6}.dim = 3;
% p_info{6}.g_type = 'cart';
% p_info{6}.name = 'SI_MIP_C=4_LS8_PWLD_XYZ=16';
% p_info{6}.title = 'Y/X = 16';
% p_info{6}.symbol = 'b--';
% 7th Info Set
% p_info{7} = struct();
% p_info{7}.dim = 3;
% p_info{7}.g_type = 'cart';
% p_info{7}.name = 'SI_MIP_C=4_LS8_PWLD_XYZ=64';
% p_info{7}.title = 'Y/X = 64';
% p_info{7}.symbol = 'k--';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin Plotting Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data
figure(1); hFig = figure(1); hold on;
set(hFig,'Position',[1,1,1200,700])
minX = 1e15; maxX = 0; minY = 1e15; maxY = 0;
for i=1:tot_sets
    p = p_info{i};
    if p.dim == 1
        d_name = [dir_root_name,num2str(p.dim),'D/'];
    else
        d_name = [dir_root_name,num2str(p.dim),'D/',p.g_type,'/'];
    end
    s_name = [d_name,p.name];
    tmat = dlmread([s_name, '.dat']);
    x = tmat(:,1); y = tmat(:,2);
    plot(x,y,p.symbol,'LineWidth',2.0)
    if min(x) < minX, minX = min(x); end
    if max(x) > maxX, maxX = max(x); end
    if min(y) < minY, minY = min(y); end
    if max(y) > maxY, maxY = max(y); end
end
% adjust legend information
for q=1:tot_sets
    legendInfo{q} = p_info{q}.title;
end
legend(legendInfo,'Location',legend_loc);
% adjust axis/font information
if strcmp(xaxis_lim, 'auto')
    xlim([minX,maxX]);
elseif strcmp(xaxis_lim, 'manual')
    xlim(xaxis_b);
end
yax = ylim; yax(1) = 0; ylim(yax);
if strcmp(yaxis_lim, 'auto')
    if maxY >= 10
        set(gca,'yscale','log');
        logminY = floor(log10(minY)); logmaxY = ceil(log10(maxY));
        ylim([10^(logminY), 10^(logmaxY)]);
    end
elseif strcmp(yaxis_lim, 'manual')
    yax = ylim; yax(2) = yaxis_max;
    ylim(yax);
end
box on;
set(gca,'xscale','log');
set(gca,'XGrid','on','XMinorGrid','off');
set(gca,'YGrid','on','YMinorGrid','off');
set(gca,'FontName','Times New Roman','FontSize',font_size);
xlabel('X in MFP', 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
ylabel('Spectral Radius', 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
% save plots
f_name = [out_name, save_name];
savefig(hFig,f_name)
export_fig(f_name,'-transparent')