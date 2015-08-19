function [X,k] = sqrtm_pulay(A,D)
%SQRTM_PULAY   Matrix square root by Pulay iteration.
%   [X,K] = SQRTM_PULAY(A,D) computes the principal square root of the
%   matrix A using the Pulay iteration with diagonal matrix D
%   (default: D = DIAG(DIAG(A))).  D must have positive diagonal entries.
%   K is the number of iterations.
%   Note: this iteration is linearly converent and converges only when
%         SQRTM(D) sufficiently well approximates SQRTM(A).

if nargin < 2, D = diag(diag(A)); end

if any(diag(D)<0) || ~isreal(D)
   error('D must have positive, real diagonal.')
end

n = length(A);
dhalf = sqrt(diag(D));
Dhalf = diag(dhalf);
B = zeros(n);
maxit = 50;

tol = mft_tolerance(A);

for k = 1:maxit

   Bold = B;
   B = (A - D - Bold^2) ./ (dhalf(:,ones(1,n)) + dhalf(:,ones(1,n))');
   X = Dhalf + B;

   reldiff = norm(B - Bold,inf)/norm(X,inf);

   if reldiff <= tol, return, end

end
error('Not converged after %2.0f iterations', maxit)
