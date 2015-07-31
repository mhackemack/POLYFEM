function out = transform_3d_plane_to_2d(v)
%assume points are Nx3, where N is number of points.
%let first point origin in 2D coordinate system (can be shifted later)
%Calculate out-of-plane vector (local z)
N = size(v,1);
origin = v(1,:);
localz = cross(v(2,:)-origin, v(3,:)-origin);
%normalize it
unitz = localz/norm(localz,2);
%calculate local x vector in plane
localx = v(2,:)-origin;   
unitx = localx/norm(localx,2);
%calculate local y
localy = cross(localz, localx);  
unity = localy/norm(localy,2); 
%assume transformation matrix
T = [localx(:), localy(:), localz(:), origin(:); 0 0 0 1];
C = [v, ones(N,1)];
out = T \ C';
out = out(1:3,:)';
return