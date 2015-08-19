function [U,H,k] = polar_newton(A)
%POLAR_NEWTON Polar decomposition by scaled Newton iteration.
%   [U,H,k] = POLAR_NEWTON(A), where the matrix A is square and
%   nonsingular, computes a unitary U and a Hermitian positive
%   definite H such that A = U*H.   k is the number of iterations.
%   A Newton iteration with acceleration parameters is used.

accel_tol = 1e-2;  % Precise value not important.
need_accel = 1;

tol = mft_tolerance(A);

X = A;
reldiff = inf;
maxit  = 16;

for k = 1:maxit

      Xold = X;
      Xinv = inv(X);
      if need_accel
         g = ( norm(Xinv,1)*norm(Xinv,inf) / (norm(X,1)*norm(X,inf)) )^(1/4);
      else
         g = 1;
      end
      X = 0.5*(g*X + Xinv'/g);
      reldiff_old = reldiff;
      diff_F = norm(X-Xold,'fro');
      reldiff = diff_F/norm(X,'fro');

      if need_accel && (reldiff < accel_tol), need_accel = false; end
      cged = (diff_F <= sqrt(tol)) || (reldiff > reldiff_old/2 && ~need_accel);
      if cged, break, end

      if k == maxit, error('Not converged after %2.0f iterations', maxit), end

end

U = X;
if nargout >= 2
   H = U'*A;
   H = (H + H')/2;  % Force Hermitian by taking nearest Hermitian matrix.
end
