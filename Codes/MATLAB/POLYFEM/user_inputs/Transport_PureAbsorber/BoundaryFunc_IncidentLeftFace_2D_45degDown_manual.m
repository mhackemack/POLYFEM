function out = BoundaryFunc_IncidentLeftFace_2D_45degDown_manual(xx,ang)
n = size(xx,1);
rang = [1,-1]/norm([1,-1]);
val = cos(pi/4);
if abs(1-dot(rang, ang/norm(ang))) < 1e-12
    out = ones(n,1)/val;
else
    out = zeros(n,1);
end