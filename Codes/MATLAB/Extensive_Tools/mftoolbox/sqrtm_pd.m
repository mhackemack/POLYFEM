function [X,k] = sqrtm_pd(A)
%SQRTM_PD    Square root of positive definite matrix via polar decomposition.
%   [X,K] = SQRT_PD(A) computes the Hermitian positive definite
%   square root X of the Hermitian positive definite matrix A.
%   It computes the Hermitian polar factor of the Cholesky factor of A.
%   K is the number of Newton polar iterations used.

R = chol(A);
[U,H,k] = polar_newton(R);
X = H;
