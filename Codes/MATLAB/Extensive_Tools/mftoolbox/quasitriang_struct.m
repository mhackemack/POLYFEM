function [m, s, k] = quasitriang_struct(R)
%QUASITRIANG_STRUCT  Block structure of upper quasitriangular matrix.
%   [M,S,K] = QUASITRIANG_STRUCT(R), where R is an upper
%   quasitriangular matrix, determines that R has M diagonal blocks,
%   the i'th of which has order S(i) and starting position K(i).
%   Any subdiagonal elements less than the tolerance EPS*NORM(R,'FRO')
%   are treated as zero.

n = length(R);
tol = eps*norm(R,'fro');

i = 1; j = 1;

while i < n
      k(j) = i;
      if abs(R(i+1,i)) <= tol
         s(j) = 1;
      else
         s(j) = 2;
      end
      i = i + s(j);
      j = j+1;
end

if i == n
   k(j) = n;
   s(j) = 1;
   m = j;
else
   m = j-1;
end
