function c = logm_cond(A)
%LOGM_COND  Relative condition number of matrix logarithm.
%   LOGM_COND(A) is the relative condition number in the Frobenius
%   norm of the matrix logarithm at the matrix A.

n = length(A);
N = n^2;
K = zeros(N);
E = zeros(n);

if ~isequal(A,triu(A))
   A = schur(A,'complex');
end

for j = 1:N
    e = zeros(N,1); e(j) = 1;
    E(:) = e;
    X = logm_frechet_pade(A,E);
    K(:,j) = X(:);
end

c = norm(K) * norm(A,'fro') / norm(logm(A),'fro');
