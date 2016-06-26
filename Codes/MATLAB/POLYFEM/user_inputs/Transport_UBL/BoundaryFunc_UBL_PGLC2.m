function out = BoundaryFunc_UBL_PGLC2(xx,ang)
n = size(xx,1);
rang1 = [3.399810435848563e-01,6.649860487269459e-01];
rang2 = [3.399810435848563e-01,-6.649860487269459e-01];
wt = 1.024387213795176e+00;
if abs(1-dot(rang1/norm(rang1), ang/norm(ang))) < 1e-12
    out = ones(n,1)/2/wt;
elseif abs(1-dot(rang2/norm(rang2), ang/norm(ang))) < 1e-12
    out = ones(n,1)/2/wt;
else
    out = zeros(n,1);
end