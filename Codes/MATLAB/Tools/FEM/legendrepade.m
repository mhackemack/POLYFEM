function [ahat,bhat]=legendrepade(f,tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% legendrepade.m
%
% Legendre-Pade Algorithm
%
% This script takes the function f(x) sampled at the Legendre-Gauss-Lobatto 
% nodes and computes the legendre coefficients for the functions a(x) and b(x)
% such that a(x)~=b(x)f(x)
%
% Note: the number of nodes (length of f) must be odd.
%
% Parameters: f    - function sampled at LGL nodes
%             ahat - expansion coeffs of numerator
%             bhat - expansion coeffs of denominator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nc=length(f);
N=(Nc-1)/2;
N1=N+1;

[x,w,P]=lglnodes(Nc-1);

% Convert to Legendre Coefficients
fhat=P\f;

% Legendre to Powers matrix
M = zeros(Nc);
M(1,1) = 1;
M(2,2) = 1;

for k=1:Nc-2
    M(:,k+2) = ((2*k+1)*[0; M(1:Nc-1,k+1)] - k*M(:,k))/(k+1);
end

% Convert to powers (exploiting parity)
fe=fhat(1:2:Nc);
fo=fhat(2:2:Nc);

Me=M(1:2:Nc,1:2:Nc);
Mo=M(2:2:Nc,2:2:Nc);

ce=Me*fe;
co=Mo*fo;

C=zeros(Nc,1);
C(1:2:Nc)=ce;
C(2:2:Nc)=co;

% Do Pade expansion from powers
lhs=eye(N1);
lhs(2:N1,2:N1)=toeplitz(C(N1:Nc-1),C(N1:-1:2));
rhs=[1;-C(N1+1:Nc)];

% Find denominator coefficients of power series
den=qmr(lhs,rhs,tol,20);

% Compute the numerator coefficients of power series
lhs=toeplitz(C(1:N1),[C(1) zeros(1,N)]);
num=lhs*den;

% Compute powers-to-Legendre matrix
X=zeros(Nc,Nc);
X(:,1)=1;
X(:,2)=x;

for k=3:Nc
    X(:,k)=x.*X(:,k-1);
end

% Point space numerator and denominator
a=X*[num;zeros(N,1)];
b=X*[den;zeros(N,1)];

% modal space numerator and denominator
ahat=P\a;
bhat=P\b;
