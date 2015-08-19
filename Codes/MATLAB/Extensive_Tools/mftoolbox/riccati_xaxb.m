function X = riccati_xaxb(A,B)
%RICCATI_XAXB  Solve Riccati equation XAX = B in positive definite matrices.
%   X = RICCATI_XAXB(A,B) is the Hermitian positive definite
%   solution to XAX = B, where A and B are Hermitian positive
%   definite matrices.

R = chol(A);
S = chol(B);

U = polar_newton(S*R');
X = R\(U'*S);
X = (X + X')/2;

% [U,H] = polar_newton(R*S'); % Variant derived in solution of Problem 6.21.
% X = R\(U*S);
