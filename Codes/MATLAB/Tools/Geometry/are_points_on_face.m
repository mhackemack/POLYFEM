%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Point on Face Calculator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to determine if a set of points reside on a
%                   geometric surface.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fbool, w_face] = are_points_on_face(verts, faces, p)
[nv, dim] = size(verts);
nf = length(faces);
np = size(p,1);
% Allocate Memory
fbool = logical(zeros(np,1));
w_face = zeros(np, 1);
% Calculate point-on-face
if dim == 1
    % Test first vertex
    fbool(abs(p - verts(1)) < 1e-13) = true;
    w_face(abs(p - verts(1)) < 1e-13) = 1;
    % Test second vertex
    fbool(abs(p - verts(2)) < 1e-13) = true;
    w_face(abs(p - verts(2)) < 1e-13) = 2;
elseif dim == 2
    [~, fbool] = inpoly(p, verts);
    for i=1:np
        if ~fbool(i), continue; end
        x = p(i,:);
        for f=1:nf
            fv = faces{f};
            L = norm(diff(verts(fv,:)));
            L1 = norm(x - verts(fv(1),:));
            L2 = norm(x - verts(fv(2),:));
            if abs(L1+L2-L) < 1e-13
                w_face(i) = f;
                break;
            end
        end
    end
elseif dim == 3
    for i=1:np
        x = p(i,:);
        phi_sum = 0;
        for f=1:nf
            fv = faces{f}; nfv = length(fv);
            for j=1:nfv
                jj = [j,mod(j,nfv)+1];
                vt = [verts(jj,:);x];
                p1 = vt(1,:) - x;
                p2 = vt(2,:) - x;
                n1 = norm(p1);
                n2 = norm(p2);
                cphi = (p1*p2') / (n1*n2);
                phi_sum = phi_sum + acos(cphi);
            end
            if abs(2*pi - phi_sum) < 1e-14
                fbool(i) = true;
                w_face(i) = f;
                break;
            end
        end
    end
end