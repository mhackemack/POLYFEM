function [S,N] = signm(A)
%SIGNM   Matrix sign decomposition.
%   [S,N] = SIGNM(A) is the matrix sign decomposition A = S*N,
%   computed via the Schur decomposition.
%   S is the matrix sign function, sign(A).

[Q, T] = schur(A,'complex');
S = Q * matsignt(T) * Q';

if nargout == 2
   N = S*A;
end

%%%%%%%%%%%%%%%%%%%%%%%%
function S = matsignt(T)
%MATSIGNT    Matrix sign function of a triangular matrix.
%   S = MATSIGN(T) computes the matrix sign function S of the
%   upper triangular matrix T using a recurrence.

n = length(T);
S = diag( sign( diag(real(T)) ) );
for p = 1:n-1
   for i = 1:n-p

      j = i+p;
      d = T(j,j) - T(i,i);

      if S(i,i) ~= -S(j,j)  % Solve via S^2 = I if we can.

         % Get S(i,j) from S^2 = I.
         k = i+1:j-1;
         S(i,j) = -S(i,k)*S(k,j) / (S(i,i)+S(j,j));

      else

         % Get S(i,j) from S*T = T*S.
         s = T(i,j)*(S(j,j)-S(i,i));
         if p > 1
            k = i+1:j-1;
            s = s + T(i,k)*S(k,j) - S(i,k)*T(k,j);
         end
         S(i,j) = s/d;

      end

   end
end
