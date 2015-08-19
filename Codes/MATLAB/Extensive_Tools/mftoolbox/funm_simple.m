function F = funm_simple(A,fun)
%FUNM_SIMPLE Simplified Schur-Parlett method for function of a matrix.
%   F = FUNM_SIMPLE(A,FUN) evaluates the function FUN at the
%   square matrix A by the Schur-Parlett method using the scalar
%   Parlett recurrence (and hence without blocking or reordering).
%   This function is intended for matrices with distinct eigenvalues
%   only and can be numerically unstable.
%   FUNM should in general be used in preference.

n = length(A);

[Q,T] = schur(A,'complex');   % Complex Schur form.
F = diag(feval(fun,diag(T))); % Diagonal of F.

% Compute off-diagonal of F by scalar Parlett recurrence.
for j=2:n
   for i = j-1:-1:1
      s = T(i,j)*(F(i,i)-F(j,j));
      if j-i >= 2
         k = i+1:j-1;
         s = s + F(i,k)*T(k,j) - T(i,k)*F(k,j);
      end
      d = T(i,i) - T(j,j);
      if d ~= 0
         F(i,j) = s/d;
      end
   end
end

F = Q*F*Q';
