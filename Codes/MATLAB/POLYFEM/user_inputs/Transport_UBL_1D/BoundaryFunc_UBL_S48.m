function out = BoundaryFunc_UBL_S48(xx,ang)
n = size(xx,1);
rang = 3.238017096286933e-02;
wt = 6.473769681268390e-02;
% val = dot(rang1(1),1);
if abs(rang - ang) < 1e-12
    out = ones(n,1)/wt;
else
    out = zeros(n,1);
end