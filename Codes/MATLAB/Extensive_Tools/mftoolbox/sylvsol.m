function X = sylvsol(A,B,C)
%SYLVSOL  Solve Sylvester equation.
%   X = SYLSOL(A, B, C) is the solution to the Sylvester equation
%   AX + XB = C, where A is m-by-m, B is n-by-n and C is m-by-n.
%   Schur decompositions are used to convert to a (quasi)-triangular
%   system.

%   Reference:
%   R. H. Bartels and G. W. Stewart.
%   Algorithm 432: Solution of the matrix equation AX+XB=C.
%   Comm. ACM, 15(9):820-826, 1972.

[m,m] = size(A);
[n,n] = size(B);

realdata = (isreal(A) && isreal(B) && isreal(C));
if ~isequal(A,triu(A)) || ~isequal(B,triu(B))

   [Q, T] = schur(A,'complex');
   [P, U] = schur(B,'complex');
   C = Q'*C*P;
   schur_red = 1;

else

   schur_red = 0;
   T = A; U = B;

end

X = zeros(m,n);

% Forward substitution.
for i = 1:n
    X(:,i) = (T + U(i,i)*eye(m)) \ (C(:,i) - X(:,1:i-1)*U(1:i-1,i));
end

if schur_red, X = Q*X*P'; end
if realdata, X = real(X); end
