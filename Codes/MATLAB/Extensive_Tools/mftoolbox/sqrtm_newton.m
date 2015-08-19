function [X,k] = sqrtm_newton(A,scal,maxit)
%SQRTM_NEWTON  Matrix square root by Newton iteration (unstable).
%   [X,k] = SQRTM_NEWTON(A,SCAL,MAXIT) computes the principal
%   square root X of the matrix A by an unstable Newton iteration.
%   SCAL specifies whether scaling is used (SCAL = 1) or not
%   (SCAL = 0, default).
%   MAXIT (default 25) is the maximum number of iterations allowed.
%   k is the number of iterations.

n = length(A);
if nargin < 2, scal = 0; end
if nargin < 3, maxit = 25; end

tol = mft_tolerance(A);
accel_tol = 1e-2;  %  Precise value not important.
need_accel = 1;

X = A;

for k = 1:maxit

   if need_accel && scal > 0
      mu = (abs(det(X))/sqrt(abs(det(A))))^(-1/n);
      X = mu*X;
   end

   Xold = X;
   X = (X + X\A)/2;

   reldiff = norm(X-Xold,'fro')/norm(X,'fro');
   if need_accel && (reldiff < accel_tol), need_accel = false; end

   if reldiff <= tol, return, end

end
error('Not converged after %2.0f iterations', maxit)
