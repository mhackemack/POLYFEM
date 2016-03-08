function out = BoundaryFunc_IncidentLeftFace_2D_45degUp_LS4(xx,ang)
n = size(xx,1);
rang = [1,1]/norm([1,1]);
val = cos(pi/4);
wt = 5.235987755982988e-01;
if abs(1-dot(rang, ang/norm(ang))) < 1e-12
    out = ones(n,1)/val/wt;
else
    out = zeros(n,1);
end