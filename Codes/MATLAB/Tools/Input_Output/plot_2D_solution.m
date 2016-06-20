%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Plot 2D Solution
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to plot scalar solutions in 2D.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_2D_solution(mesh, DoF, ~, x, basis_name)
n_input = nargin;
if iscell(x), x = x{1}; end
% determine plotting type
if n_input < 5
    if DoF.DoFType == 0
        plot_LD_solution(mesh, DoF, x);
    elseif DoF.Degree == 0
        plot_k0_solution(mesh,DoF,x);
    else
        plot_general_solution(mesh, DoF, x);
    end
else
    basis_name = upper(basis_name);
    if strcmp(basis_name, 'PWLD')
        plot_general_PWLD_solution(mesh, DoF, x);
    elseif strcmp(basis_name, 'PWQ')
        
    else
        plot_general_solution(mesh, DoF, x);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_LD_solution(mesh,DoF,x)
hold on
for c=1:mesh.TotalCells
    cv = mesh.CellVerts{c};
    xx = mesh.Vertices(cv,:);
    ccent = mesh.CellCenter(c,:);
    dx = get_LD_widths(mesh.Dimension,xx);
    b = get_LD_basis(xx,ccent,dx);
    cn = DoF.ConnectivityArray{c};
    y    = b*x(cn);
    patch(xx(:,1), xx(:,2), y, y, 'FaceColor', 'interp')
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_k0_solution(mesh,DoF,x)
hold on
for c=1:mesh.TotalCells
    cv = mesh.CellVerts{c}; ncv = length(cv);
    xx = mesh.Vertices(cv,:);
    cn = DoF.ConnectivityArray{c};
    y  = x(cn)*ones(ncv,1);
    patch(xx(:,1), xx(:,2), y, y, 'FaceColor', 'interp')
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_general_solution(mesh,DoF,x)
hold on
for c=1:mesh.TotalCells
    cv = mesh.CellVerts{c}; ncv = length(cv);
    n = DoF.ConnectivityArray{c}; n(ncv+1:end) = [];
    nlocs = DoF.NodeLocations(n,:);
    patch(nlocs(:,1), nlocs(:,2), x(n), x(n), 'FaceColor', 'interp')
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_general_PWLD_solution(mesh, DoF, x)
hold on
for c=1:mesh.TotalCells
    cv = mesh.get_cell_verts(c);
    xc = mean(cv);
    ncv = size(cv,1);
    n = DoF.ConnectivityArray{c}; n(ncv+1:end) = [];
    xx = x(n);
    xxc = mean(xx);
    for j=1:ncv
        if j==ncv
            jj = [j,1];
        else
            jj = [j,j+1];
        end
        xxx = [xx(jj);xxc];
        patch([cv(jj,1);xc(1)], [cv(jj,2);xc(2)], xxx, xxx, 'FaceColor', 'interp', 'LineStyle', 'none')
    end
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xmean,dx] = get_LD_widths(dim,verts)
xmean = mean(verts);
% xmean = compute_volume_centroid(dim, verts, faces);
if dim == 1
    dx = verts(2) - verts(1);
elseif dim == 2
    dx = [max(verts(:,1)) - min(verts(:,1)),...
          max(verts(:,2)) - min(verts(:,2))];
elseif dim == 3
    dx = [max(verts(:,1)) - min(verts(:,1)),...
          max(verts(:,2)) - min(verts(:,2)),...
          max(verts(:,3)) - min(verts(:,3))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bout = get_LD_basis(x,xmean,dx)
nx = size(x,1); onx = ones(nx,1);
ddr = xmean./dx;
bout = [onx, 2*x./(onx*dx) - 2*onx*ddr];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%