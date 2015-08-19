function [X,k] = signm_newton(A,scal)
%SIGNM_NEWTON   Matrix sign function by Newton iteration.
%   [S,k] = SIGNM_NEWTON(A,SCAL) computes the matrix sign
%   function S of A using the scaled Newton iteration, with
%   scale factors specified by SCAL:
%   SCAL = 0: no scaling.
%   SCAL = 1: determinantal scaling (default).
%   SCAL = 2: spectral scaling.
%   SCAL = 3: norm scaling.
%   k is the number of iterations.

if nargin < 2, scal = 1; end

tol = mft_tolerance(A);
accel_tol = 1e-2;  %  Precise value not important.
need_accel = 1;
n = length(A);
maxit = 16;

X = A;
reldiff = inf;

for k = 1:maxit

      Xold = X;
      Xinv = inv(X);
      if need_accel && scal > 0
         % In practice should estimate spectral radius and 2-norms;
         % here they are computed exactly.
         switch scal
         case 1
              g = abs(det(X))^(-1/n);
         case 2
              s1 = max(abs(eig(Xinv)));
              s2 = max(abs(eig(X)));
              g = sqrt(s1/s2);
         case 3
              g = sqrt( norm(Xinv) / (norm(X)) );
         end
         X = g*X; Xinv = Xinv/g;
      end
      X = 0.5*(X + Xinv);
      diff_F = norm(X-Xold,'fro');
      reldiff_old = reldiff;
      reldiff = diff_F/norm(X,'fro');

      if need_accel && (reldiff < accel_tol), need_accel = false; end
      cged = (diff_F <= sqrt( tol*norm(X)/norm(Xinv) ) || ...
              reldiff > reldiff_old/2 && ~need_accel);
      if cged, return, end

end
error('Not converged after %2.0f iterations', maxit)
