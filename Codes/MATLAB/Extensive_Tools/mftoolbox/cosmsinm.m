function [C,S,cost,s] = cosmsinm(A)
%COSMSINM  Matrix cosine and sine by double angle algorithm.
%   [C,S,COST,s] = COSMSINM(A) computes the cosine C and the sine S
%   of the matrix A using the double angle algorithm with Pade
%   approximation.  The total number of matrix multiplications and
%   linear systems solved is returned as COST and s specifies the
%   amount of scaling (A is scaled to 2^(-s)*A).

n = length(A);
theta = [0.00613443965526 0.11110098037055  0.43324697422211 0.98367255274344 ...
         1.72463663220280 2.61357494421368 3.61521023400301 4.70271938553349 ...
         5.85623410320942 7.06109053959248 0 9.58399601173102];

normA = norm(A,inf);
d = [2 4 6 8 10 12 14 16 20];
for i = d(1:8)
    if normA <= theta(i/2);
        m = i;
        cost = m/2+3;
        [C,S] = cosmsinm_pade(A,m);
        s = m;
        return
    end
end

s = 0;
if normA > theta(20/2)
    s = max( ceil( log2( normA/theta(20/2) ) ), 1 );
    A = A/(2^s);
    normA = normA/(2^s);
end

if normA > 2*theta(12/2)
    m = 20; cost = 12;
elseif normA > theta(16/2)
    A = A/2; m = 12; cost = 9; s = s+1;
else
    m = 16; cost = 11;
end

[C,S] = cosmsinm_pade(A,m);
for i = 1:s
    S = 2*S*C;
    C = 2*(C^2) - eye(n);
end
cost = cost + 2*s;
