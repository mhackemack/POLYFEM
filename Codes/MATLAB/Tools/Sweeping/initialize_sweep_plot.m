%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Initialize Sweep Plot
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_sweep_plot(mesh)
dim = mesh.Dimension;
if dim == 1
    plot_mesh(mesh,0,0);
elseif dim == 2
    plot_mesh(mesh,0,0);
elseif dim == 3
    plot_mesh(mesh,0,0,0);
    view(-35,45);
end
hold on;