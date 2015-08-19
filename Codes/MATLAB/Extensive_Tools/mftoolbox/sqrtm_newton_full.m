function [X,k] = sqrtm_newton_full(A, X0)
%SQRTM_NEWTON_FULL   Matrix square root by full Newton method.
%   [X,K] = SQRTM_NEWTON_FULL(A, X0) applies Newton's method to
%   compute a square root X of the matrix A, with starting matrix X0.
%   Default: X0 = A.  K is the number of iterations.

if nargin < 2, X0 = A; end

X = X0;
tol = mft_tolerance(A);
maxit = 50;

for k = 1:maxit

   Xold = X;
   R = A - X^2;
   % Solve XE + EX = R.
   E = sylvsol(X,X,R);
   X = X + E;

   reldiff = norm(X - Xold,inf)/norm(X,inf);

   if reldiff <= tol; return; end

end
error('Not converged after %2.0f iterations', maxit)
