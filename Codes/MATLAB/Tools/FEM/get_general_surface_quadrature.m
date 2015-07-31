%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          General Quadrature Generator - Surface
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce a quadrature set for general
%                   polygons/polyhedra.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_general_surface_quadrature(verts, f_ind, deg, poly_bool)
% Get General Input Information 
% -----------------------------
if nargin < 4, poly_bool = false; end
[nverts, dim] = size(verts);
fverts = verts(f_ind,:);
nf_ind = length(f_ind);
% Quick 1D Output
% ---------------
if dim == 1
    qw = 1; qx = verts(f_ind);
    return
end
% Calculate Quadrature Set
% ------------------------
[rqx, rqw] = get_reference_quadrature(dim, nf_ind, deg); nrqx = length(rqw);
if dim == 2
    v0 = fverts(1,:); df = diff(verts(f_ind,:));
    flen = norm(df);
    qw = rqw * flen;
    qx = [df(1)*rqx+v0(1), df(2)*rqx+v0(2)];
else
    farea = polygonArea3d(fverts);
    if nf_ind == 3
        qw = rqw / sum(rqw) * farea;
        
    elseif nf_ind == 4 % isoparametric mapping on the 3d quad (not tested for concave faces)
        qx = ((1-rqx(:,1)).*(1-rqx(:,2)))*fverts(1,:) + (rqx(:,1).*(1-rqx(:,2)))*fverts(2,:)...
             + (rqx(:,1).*rqx(:,2))*fverts(3,:) + (rqx(:,2).*(1-rqx(:,1)))*fverts(4,:);
        qw = rqw / sum(rqw) * farea;
    else
        qx = zeros(nrqx*nf_ind,dim); qw = zeros(nrqx*nf_ind,1);
        fmean = mean(verts(f_ind,:));
        for i=1:nf_ind
            if i==nf_ind
                ii = [i,1];
            else
                ii = [i,i+1];
            end
            tv = [fverts(ii,:);fmean];
            ind = (i-1)*nrqx+1:i*nrqx;
            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_reference_quadrature(dim, nf_ind, deg)
if dim == 2
    [qx, qw] = get_legendre_gauss_quad(deg);
else
    if nf_ind == 4
        [tqx, tqw] = get_legendre_gauss_quad(deg); nt = length(tqw);
        qw = tqw*tqw'; qw = qw(:); qw = qw / sum(qw);
        qxx = zeros(nt, nt); qxy = zeros(nt, nt);
        for i=1:nt
            qxx(i,:) = tqx(i); qxy(:,i) = tqx(i);
        end
        qx = [qxx(:), qxy(:)];
    else
        [qx, qw] = Quad_On_Triangle(deg);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%