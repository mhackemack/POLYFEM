function [x,w,C]=gglnodes(N,a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% gglnodes.m
%
% Computes the Gegenbauer-Gauss-Lobatto nodes, weights and the GGL 
% Vandermonde matrix. The GGL nodes are the zeros of (1-x^2)*C'_N(x). 
% Useful for numerical integration and spectral methods. 
%
% Written by Greg von Winckel - 05/27/2005
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation + 1
N1=N+1;

x=cos(pi*(0:N)'/N);

% The Gegenbauer Vandermonde Matrix
C=zeros(N1,N1);

% Compute G_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;

if a>1/2

asteps=ceil(2*a);
da=a/asteps;
   
    
for step=1:asteps
    
A=da*step;
    
    while max(abs(x-xold))>eps

        xold=x;
    
        C(:,1)=1;    C(:,2)=2*A*x;
    
        for k=2:N
            C(:,k+1)=(2*(k-1+A)*x.*C(:,k)-(k+2*A-2)*C(:,k-1))/k;
        end

        num=-N*x.*C(:,N1)+(N+2*A-1)*C(:,N);
        den=N*(N*(x.^2-1)-2*A+1).*C(:,N1)+(2*N*A-N+4*A*(A-1))*x.*C(:,N);

        x=xold-(1-x.^2).*num./den;

        end
        xold=2;

    end

else
    
    while max(abs(x-xold))>eps
    
        xold=x;
    
        C(:,1)=1;    C(:,2)=2*a*x;
    
        for k=2:N
            C(:,k+1)=(2*(k-1+a)*x.*C(:,k)-(k+2*a-2)*C(:,k-1))/k;
        end

        num=-N*x.*C(:,N1)+(N+2*a-1)*C(:,N);
        den=N*(N*(x.^2-1)-2*a+1).*C(:,N1)+(2*N*a-N+4*a*(a-1))*x.*C(:,N);

        x=xold-(1-x.^2).*num./den;
    end
end


% Compute weights
w=pi^(-1/2)*(C')\[gamma(1/2+a)/gamma(1+a);zeros(N,1)];
