function out = BoundaryFunc_UBL_S24(xx,ang)
n = size(xx,1);
rang = 6.405689286260557e-02;
wt = 1.279381953467523e-01;
% val = dot(rang1(1),1);
if abs(rang - ang) < 1e-12
    out = ones(n,1)/wt;
else
    out = zeros(n,1);
end