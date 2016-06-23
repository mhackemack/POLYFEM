function out = BoundaryFunc_UBL_PGLC8(xx,ang)
n = size(xx,1);
rang1 = [9.501250983763748e-02,7.039078856549176e-01];
rang2 = [9.501250983763748e-02,-7.039078856549176e-01];
wt = 2.975883230118723e-01;
% val = dot(rang1(1),1);
if abs(1-dot(rang1/norm(rang1), ang/norm(ang))) < 1e-12
    out = ones(n,1)/2/wt;
elseif abs(1-dot(rang2/norm(rang2), ang/norm(ang))) < 1e-12
    out = ones(n,1)/2/wt;
else
    out = zeros(n,1);
end