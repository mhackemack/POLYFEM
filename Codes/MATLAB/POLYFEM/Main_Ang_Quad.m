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
clear; clc; close all; format long e;
out_dir = 'outputs/Ang_Quads/';

alvls = [2,4,6,8];
plvls = [4,6,8];
snlvl = [2,4,8,16,24];

dim = 3;
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.QuadType = 'PGLC';
fname = data.Neutronics.Transport.QuadType;

if strcmp(data.Neutronics.Transport.QuadType, 'PGLC')
    for i=1:length(plvls)
        data.Neutronics.Transport.PolarLevels = plvls(i);
        for j=1:length(alvls)
            data.Neutronics.Transport.AzimuthalLevels = alvls(j);
            data.Neutronics.Transport = get_angular_quadrature(data.Neutronics.Transport, dim);
            draw_quadrature(data);
            ffname = [out_dir,fname,num2str(alvls(j)),'_',num2str(plvls(i))];
%             export_fig ffname -eps
            saveas(gcf,ffname,'png');
            close all;
        end
    end
else
    for m=1:length(snlvl)
        data.Neutronics.Transport.SnLevels = snlvl(m);
        data.Neutronics.Transport = get_angular_quadrature(data.Neutronics.Transport, dim);
        draw_quadrature(data);
    end
end