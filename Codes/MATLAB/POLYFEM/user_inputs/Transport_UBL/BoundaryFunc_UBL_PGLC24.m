function out = BoundaryFunc_UBL_PGLC24(xx,ang)
n = size(xx,1);
rang1 = [3.238017096286933e-02,7.067359919122612e-01];
rang2 = [3.238017096286933e-02,-7.067359919122612e-01];
% wt = (1.016897363585255e-01)/2;
% val = dot(rang1(1),1);
if abs(1-dot(rang1/norm(rang1), ang/norm(ang))) < 1e-12
    out = ones(n,1);
elseif abs(1-dot(rang2/norm(rang2), ang/norm(ang))) < 1e-12
    out = ones(n,1);
else
    out = zeros(n,1);
end