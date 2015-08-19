function [P,Q,k] = sqrtm_db(A,scale)
%SQRTM_DB   Matrix square root by Denman-Beavers iteration.
%   [P,Q,k] = SQRTM_DB(A,SCAL) computes the principal square root
%   P of the matrix A using the Denman-Beavers iteration.
%   It also returns Q = INV(P).
%   SCAL specifies the scaling:
%   SCAL == 0, no scaling.
%   SCAL == 1, determinantal scaling (default).
%   k is the number of iterations.

n = length(A);
if nargin < 2, scale = 1; end

tol = mft_tolerance(A);
P = A;
Q = eye(n);
reldiff = inf;
maxit = 25;

for k = 1:maxit

   if scale == 1
       g = (abs(det(P)*det(Q)))^(-1/(2*n));
       P = g*P; Q = g*Q;
   end

   Pold = P;

   Poldinv = inv(Pold);
   P = (P + inv(Q))/2;
   Q = (Q + inv(Pold))/2;

   diff_F = norm(P-Pold,'fro');
   reldiff_old = reldiff;
   reldiff = diff_F/norm(P,'fro');
   if reldiff < 1e-2, scale = 0; end  % Switch to no scaling.

   cged = (diff_F <= sqrt( tol*norm(P)/norm(Poldinv) ) || ...
           reldiff > reldiff_old/2 && ~scale);
   if cged, return, end

end
error('Not converged after %2.0f iterations', maxit)
