function [c,K] = expm_cond(A)
%EXPM_COND  Relative condition number of matrix exponential.
%   EXPM_COND(A) is the relative condition number in the Frobenius
%   norm of the matrix exponential at the matrix A.
%   [C,K] = EXPM_COND(A) returns the condition number C and the Kronecker
%   matrix form K of the Frechet derivative.

n = length(A);
N = n^2;
K = zeros(N);
E = zeros(n);

if nargout < 2 && ~isequal(A,triu(A))
   % If returning K cannot use Schur form.
   A = schur(A,'complex');
end

for j = 1:N
    e = zeros(N,1); e(j) = 1;
    E(:) = e;
    X = expm_frechet_pade(A,E);
    K(:,j) = X(:);
end

c = norm(K) * norm(A,'fro') / norm(expm(A),'fro');
