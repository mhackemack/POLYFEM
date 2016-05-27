%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Generate Angular Quadrature Set Plots
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
% -------------------
if exist('pbool', 'var')
    clearvars -except pbool
else
    clear; pbool = false;
end
clc; close all; format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end

alvls = [1];
plvls = [12];
snlvl = [2,4,8,16];
pdim  = 1;

dim = 2;
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.QuadType = 'PGLC';
fname = data.Neutronics.Transport.QuadType;
out_dir = 'outputs/Ang_Quads/';

if strcmp(data.Neutronics.Transport.QuadType, 'PGLC')
    for i=1:length(plvls)
        data.Neutronics.Transport.PolarLevels = plvls(i);
        for j=1:length(alvls)
            data.Neutronics.Transport.AzimuthalLevels = alvls(j);
            data.Neutronics.Transport.PolarDimension = pdim;
            data.Neutronics.Transport = get_angular_quadrature(data.Neutronics.Transport, dim);
            draw_quadrature(data);
            ffname = [out_dir,fname,num2str(alvls(j)),'_',num2str(plvls(i)),'_',num2str(dim),'D'];
            saveas(gcf,ffname,'fig');
            print(gcf,ffname,'-depsc');
            print(gcf,ffname,'-dpng');
            close all;
        end
    end
else
    for m=1:length(snlvl)
        data.Neutronics.Transport.SnLevels = snlvl(m);
        data.Neutronics.Transport = get_angular_quadrature(data.Neutronics.Transport, dim);
        draw_quadrature(data);
        ffname = [out_dir,fname,num2str(snlvl(m)),'_',num2str(dim),'D'];
        saveas(gcf,ffname,'fig');
        print(gcf,ffname,'-depsc');
        print(gcf,ffname,'-dpng');
        close all;
    end
end