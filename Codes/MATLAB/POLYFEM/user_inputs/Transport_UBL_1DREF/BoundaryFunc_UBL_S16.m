function out = BoundaryFunc_UBL_S16(xx,ang)
n = size(xx,1);
rang = 9.501250983763748e-02;
wt = 1.894506104550684e-01;
% val = dot(rang1(1),1);
if abs(rang - ang) < 1e-12
    out = ones(n,1)/wt;
else
    out = zeros(n,1);
end