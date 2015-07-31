function [ahat,bhat]=chebyshevpade(f,tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% chebyshevpade.m
%
% Chebyshev-Pade Algorithm
%
% This script takes the function f(x) sampled at the Chebyshev-Gauss-Lobatto 
% nodes and computes the legendre coefficients for the functions a(x) and b(x)
% such that a(x)~=b(x)f(x)
%
% Note: the number of nodes (length of f) must be odd.
%
% Parameters: f    - function sampled at CGL nodes
%             ahat - expansion coeffs of numerator
%             bhat - expansion coeffs of denominator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nc=length(f);
Nc1=Nc-1;
N=(Nc-1)/2;
N1=N+1;

x=cos(pi*(0:Nc1)/(Nc1))';

F=real(ifft([f(1:Nc);f(Nc1:-1:2)]));
fhat=[F(1,:); 2*F(2:Nc,:)];
    
% Legendre to Powers matrix
% Chebyshev to Powers Matrix
M = zeros(Nc);
M(1,1) = 1;
M(2,2) = 1;

for k=1:Nc1-1
    M(:,k+2) = 2*[0; M(1:Nc1,k+1)] - M(:,k);
end

% Convert to powers (exploiting parity)
fe=fhat(1:2:Nc);
fo=fhat(2:2:Nc);

Me=M(1:2:Nc,1:2:Nc);    Mo=M(2:2:Nc,2:2:Nc);
ce=Me*fe;               co=Mo*fo;

C=zeros(Nc,1);
C(1:2:Nc)=ce;           C(2:2:Nc)=co;

% Do Pade expansion from powers
lhs=eye(N1);
lhs(2:N1,2:N1)=toeplitz(C(N1:Nc-1),C(N1:-1:2));
rhs=[1;-C(N1+1:Nc)];

% Find denominator coefficients of power series
den=qmr(lhs,rhs,tol,20);

% Compute the numerator coefficients of power series
lhs=toeplitz(C(1:N1),[C(1) zeros(1,N)]);
num=lhs*den;

% Compute powers-to-Chebyshev matrix
X=zeros(Nc,Nc);
X(:,1)=1;
X(:,2)=x;

for k=3:Nc
    X(:,k)=x.*X(:,k-1);
end

% Point space numerator and denominator
a=X*[num;zeros(N,1)];
b=X*[den;zeros(N,1)];

plot(x,f,x,a./b,'o');

% modal space numerator and denominator
temp=real(ifft([a(1:Nc);a(Nc1:-1:2)]));
ahat=[temp(1,:); 2*temp(2:Nc,:)];

temp=real(ifft([b(1:Nc);b(Nc1:-1:2)]));
bhat=[temp(1,:); 2*temp(2:Nc,:)];