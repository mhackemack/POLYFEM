function c = fab_arnoldi(A,b,fun,m)
%FAB_ARNOLDI  f(A)*b approximated by Arnoldi method.
%   C = FAB_ARNOLDI(A,B,FUN,M) approximates FUNM(A,FUN)*B
%   for a square matrix A using M steps of the Arnoldi process
%   with starting vector B/norm(B).
%   FUN must be a function handle for which FUNM(A,FUN) is defined.
%   For large matrices M is intended to be much less than LENGTH(A).

q1 = b/norm(b);
[Q,H] = arnoldi(A,q1,m);
H = H(1:m,1:m);
Q = Q(:,1:m);
e = zeros(m,1); e(1) = 1;
c = norm(b)*Q*funm(H,fun)*e;
