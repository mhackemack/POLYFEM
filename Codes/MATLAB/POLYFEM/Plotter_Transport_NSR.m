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
dir_root_name = 'outputs/Transport_NSR/';
out_name = [dir_root_name,'figures/'];
% Begin User Input Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_sets = 3; p_info = cell(tot_sets, 1);
save_name = 'MIP_V_hex_LS2,4,8_C=4';
legend_loc = 'northwest'; font_size = 22;
% 1st Info Set
p_info{1} = struct();
p_info{1}.title = 'MIP, S2';
p_info{1}.symbol = 'ko-';
p_info{1}.diff_type = 'MIP';
p_info{1}.m_type = 'hex';
p_info{1}.bc_type = 'Vacuum';
p_info{1}.q_type = 'LS';
p_info{1}.Sn = 2;
p_info{1}.C = 4;
% 2nd Info Set
p_info{2} = struct();
p_info{2}.title = 'MIP, S4';
p_info{2}.symbol = 'bs-';
p_info{2}.diff_type = 'MIP';
p_info{2}.m_type = 'hex';
p_info{2}.bc_type = 'Vacuum';
p_info{2}.q_type = 'LS';
p_info{2}.Sn = 4;
p_info{2}.C = 4;
% 3rd Info Set
p_info{3} = struct();
p_info{3}.title = 'MIP, S8';
p_info{3}.symbol = 'rv-';
p_info{3}.diff_type = 'MIP';
p_info{3}.m_type = 'hex';
p_info{3}.bc_type = 'Vacuum';
p_info{3}.q_type = 'LS';
p_info{3}.Sn = 8;
p_info{3}.C = 4;
% 4th Info Set
% p_info{4} = struct();
% p_info{4}.title = 'MIP, hexahedra';
% p_info{4}.symbol = 'ko--';
% p_info{4}.diff_type = 'MIP';
% p_info{4}.m_type = 'hex';
% p_info{4}.bc_type = 'Vacuum';
% p_info{4}.q_type = 'LS';
% p_info{4}.Sn = 8;
% p_info{4}.C = 4;
% 5th Info Set
% p_info{5} = struct();
% p_info{5}.title = 'MIP, triangular prisms';
% p_info{5}.symbol = 'bs--';
% p_info{5}.diff_type = 'MIP';
% p_info{5}.m_type = 'tri_prism';
% p_info{5}.bc_type = 'Vacuum';
% p_info{5}.q_type = 'LS';
% p_info{5}.Sn = 8;
% p_info{5}.C = 4;
% 6th Info Set
% p_info{6} = struct();
% p_info{6}.title = 'MIP, polygonal prisms';
% p_info{6}.symbol = 'rv--';
% p_info{6}.diff_type = 'MIP';
% p_info{6}.m_type = 'poly3D';
% p_info{6}.bc_type = 'Vacuum';
% p_info{6}.q_type = 'LS';
% p_info{6}.Sn = 8;
% p_info{6}.C = 4;
% End User Input Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin Plotting Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data
figure(1); hFig = figure(1); hold on;
set(hFig,'Position',[1,1,1200,700])
minX = 1e15; maxX = 0; minY = 1e15; maxY = 0;
for i=1:tot_sets
    p = p_info{i};
    dir_name = [dir_root_name, p.diff_type, '_', p.bc_type, '_', p.m_type,'/'];
    s_name = [dir_name, p.q_type, num2str(p.Sn), '_C=',num2str(p.C),'.dat'];
    tmat = dlmread(s_name);
    x = tmat(:,1); y = tmat(:,end); ynorm = tmat(:,end-1);
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
xlim([minX,maxX]);
yax = ylim; yax(1) = 0; ylim(yax);
if maxY >= 10
    set(gca,'yscale','log');
    logminY = floor(log10(minY)); logmaxY = ceil(log10(maxY));
    ylim([10^(logminY), 10^(logmaxY)]);
end
box on;
set(gca,'xscale','log');
set(gca,'XGrid','on','XMinorGrid','off');
set(gca,'YGrid','on','YMinorGrid','off');
set(gca,'FontName','Times New Roman','FontSize',font_size);
xlabel('Cell Size in MFP', 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
ylabel('Spectral Radius', 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
% save plots
f_name = [out_name, save_name];
savefig(hFig,f_name)
export_fig(f_name,'-transparent')

