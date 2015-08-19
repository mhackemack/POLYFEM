function [C,cost,s] = cosm(A)
%COSM    Matrix cosine by double angle algorithm.
%   [C,COST,s] = COSM(A) computes the cosine of the matrix A using the
%   double angle algorithm with Pade approximation.  The total number of
%   matrix multiplications and linear systems solved is returned as
%   COST and s specifies the amount of scaling (A is scaled to 2^(-s)*A).

n = length(A);
% theta_m for m=2:2:20.
theta = [0.00613443965526 0.11110098037055  0.43324697422211 0.98367255274344 ...
           1.72463663220280 2.61357494421368 3.61521023400301 4.70271938553349 ...
           5.85623410320942 7.06109053959248];
s = 0;

B = A^2;
normA2 = sqrt(norm(B,inf));
d = [2 4 6 8 12 16 20];
for i = d(1:6)
    if normA2 < theta(i/2)
        m = i;
        cost = (m<=8)*(m/2+1) + (m==12)*6;
        C = cosm_pade(B,m);
        s = m;
        return
    end
end

if normA2 > theta(20/2)
    s = ceil( log2( normA2/theta(20/2) ) );
    B = B/(4^s);
    normA2 = normA2/(2^s);
end

if normA2 > 2*theta(12/2)
    m = 20; cost = 8;
elseif normA2 > theta(16/2)
    B = B/4; m = 12; cost = 6; s = s+1;
else
    m = 16; cost = 7;
end

C = cosm_pade(B,m);
for i = 1:s
    C = 2*(C^2) - eye(n);
end
cost = cost+s;
