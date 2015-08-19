function [X,k] = rootpm_newton(A,p,c)
%ROOTPM_NEWTON  Coupled Newton iteration for matrix pth root.
%   [X,k] = ROOTPM_NEWTON(A,P,C) computes the principal
%   Pth root of the matrix A.
%   C (default 1) is a convergence parameter.
%   k is the number of iterations.

if nargin < 3, c = 1; end

n = length(A);
M = A/c^p;
X = c*eye(n);
tol = mft_tolerance(A);
maxit = 20;

relres = inf;

for k=1:maxit

   X = ( ((p+1)*eye(n) - M)/p )\X;
   M = power_binary( ((p+1)*eye(n) - M)/p, p) * M;

   relres_old = relres;
   relres = norm(M-eye(n),inf);

   if relres <= tol || relres > relres_old/2, return, end

end
error('Not converged after %2.0f iterations', maxit)
