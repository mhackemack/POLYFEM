function [X,M,k] = sqrtm_dbp(A,scale)
%SQRTM_DBP  Matrix square root by product form of Denman-Beavers iteration.
%   [X,M,k] = SQRTM_DBP(A,SCAL) computes the principal square root X
%   of the matrix A using the product form of the Denman-Beavers
%   iteration. The matrix M tends to EYE.
%   SCAL specifies the scaling:
%        SCAL == 0, no scaling.
%        SCAL == 1, determinantal scaling (default).
%   k is the number of iterations.

n = length(A);
if nargin < 2, scale = 1; end

tol = mft_tolerance(A);
X = A;
M = A;
maxit = 25;

for k = 1:maxit

   if scale == 1
       g = (abs(det(M)))^(-1/(2*n));
       X = g*X; M = g^2*M;
   end

   Xold = X; invM = inv(M);

   X = X*(eye(n) + invM)/2;
   M = 0.5*(eye(n) + (M + invM)/2);

   Mres = norm(M - eye(n),'fro');

   reldiff = norm(X - Xold,'fro')/norm(X,'fro');
   if reldiff < 1e-2, scale = 0; end  % Switch to no scaling.

   if Mres <= tol, return; end

end
error('Not converged after %2.0f iterations', maxit)
