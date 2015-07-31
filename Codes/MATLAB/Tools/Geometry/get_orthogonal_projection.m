%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get Orthogonal Projection 
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:  1) verts - Cell Vertices
%           2) faces - Face Vertices
%
%           Optional:
%           3) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_orthogonal_projection(verts, faces, varargin)
% Process Input and Output Space
[nv, dim] = size(verts);
nf = length(faces);
out = zeros(nf, 1);
if dim == 1
    out(:) = abs(verts(2) - verts(1));
    return
end
% Compute Reference Quadrature
% ----------------------------
qorder = 3;
if nargin>2, qorder=varargin{1}; end
if dim == 2
    [qx, qw] = get_legendre_gauss_quad(qorder); nqx = qorder;
elseif dim == 3
    [qx, qw] = Quad_On_Triangle(qorder); nqx = length(qw);
end
cv = get_cell_volume(dim, verts, faces, nv, nf);
% Loop through faces and calculate orthogonal projections
% -------------------------------------------------------
for f=1:nf
    fv = faces{f}; nfv = length(fv);
    fc = mean(verts(fv,:));
    if dim == 1
        out(f) = cv;
    elseif dim == 2
        v0 = verts(fv(1),:);
        vv = diff(verts(fv,:));
        len = norm(vv);
        if nv == 3
            out(f) = 2*cv/len;
        elseif nv == 4
            out(f) = cv/len;
        else
            n = [-vv(2),vv(1)];
            for q=1:nqx
                vi = v0 + qx(q)*vv;
                h = get_intersections_2D(vi, n, verts, faces);
                h = abs(max(h));
                out(f) = out(f) + qw(q)*h;
            end
        end
    elseif dim == 3
        if nf == 4
            
        elseif nf == 6 && nv == 8
            
        else
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function Listing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h, i_bool] = get_intersections_2D(v, n, verts, f)
nf = length(f);
h = zeros(nf, 1);
i_bool = logical(zeros(nf, 1));
for ff=1:nf
    fv = f{ff};
    vv = [verts(fv(1),:), verts(fv(2),:)];
    inter = intersectLineEdge([v,n],vv);
    if ~isnan(inter(1)) && ~isinf(inter(1))
        i_bool(ff) = true;
        h(ff) = norm(inter - v);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h, i_bool] = get_intersections_3D(v, n, verts, f)
nf = length(f);
h = zeros(nf, 1);
i_bool = logical(zeros(nf, 1));
for ff=1:nf
    fv = f{ff};
    vv = verts(fv,:);
    [inter, inside] = intersectLinePolygon3d([v,n], vv);
    if inside
        i_bool(ff) = true;
        h(ff) = norm(inter - v);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_cell_volume(dim, v, f, nv, nf)
if dim == 1
    out = v(2) - v(1);
elseif dim == 2
    out = polygonArea(v);
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

