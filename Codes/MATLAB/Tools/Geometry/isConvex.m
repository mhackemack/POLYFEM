%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Convex/Concave Geometry Calculator
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
function out = isConvex(verts, faces)
[nv, dim] = size(verts);
if dim == 1
    out = true;
elseif dim == 2
    out = true;
    atemp = 0;
    atot = (nv-2)*pi;
    for i=1:nv
        ii = [i,mod(i,nv)+1];
        if i==1
            ii = [nv,ii];
        else
            ii = [i-1,ii];
        end
        l12 = norm(diff(verts(ii(1:2),:)));
        l23 = norm(diff(verts(ii(2:3),:)));
        l13 = norm(diff(verts(ii([1,3]),:)));
        cphi = (l12^2 + l23^2 - l13^2) / (2*l12*l23);
        if abs(cphi + 1) < 1e-12
            phi = pi;
        else
            phi = acos(cphi);
        end
        atemp = atemp + phi;
    end
    if abs(atemp - atot) > 1e-12, out = false; end
elseif dim == 3
    
end
