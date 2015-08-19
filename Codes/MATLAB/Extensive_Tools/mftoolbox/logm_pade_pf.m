function S = logm_pade_pf(A,m)
%LOGM_PADE_PF   Evaluate Pade approximant to matrix log by partial fractions.
%   Y = LOGM_PADE_PF(A,M) evaluates the [M/M] Pade approximation to
%   LOG(EYE(SIZE(A))+A) using a partial fraction expansion.

[nodes,wts] = gauss_legendre(m);
% Convert from [-1,1] to [0,1].
nodes = (nodes + 1)/2;
wts = wts/2;

n = length(A);
S = zeros(n);

for j=1:m
    S = S + wts(j)*(A/(eye(n) + nodes(j)*A));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = gauss_legendre(n)
%GAUSS_LEGENDRE  Nodes and weights for Gauss-Legendre quadrature.
%   [X,W] = GAUSS_LEGENDRE(N) computes the nodes X and weights W
%   for N-point Gauss-Legendre quadrature.

% Reference:
% G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature
% rules, Math. Comp., 23(106):221-230, 1969.

i = 1:n-1;
v = i./sqrt((2*i).^2-1);
[V,D] = eig( diag(v,-1)+diag(v,1) );
x = diag(D);
w = 2*(V(1,:)'.^2);
