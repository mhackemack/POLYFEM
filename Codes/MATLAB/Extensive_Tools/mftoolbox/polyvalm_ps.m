function [P,s,cost] = polyvalm_ps(c,A,s)
%POLYVALM_PS  Evaluate polynomial at matrix argument by Paterson-Stockmeyer alg.
%   [P,S,COST] = POLYVALM_PS(C,A,S) evaluates the polynomial whose
%   coefficients are the vector C at the matrix A using the
%   Paterson-Stockmeyer algorithm.  If omitted, the integer parameter
%   S is chosen automatically and its value is returned as an
%   output argument.   COST is the number of matrix multiplications used.

m = length(c)-1; % Degree of poly.
c = c(end:-1:1); c = c(:);
n = length(A);

if nargin < 3
   % Determine optimum parameter s.
   s = ceil(sqrt(m));
end
r = floor(m/s);
cost = s+r-(m==r*s)-1;

% Apower{i+1} = A^i;
Apower = cell(s+1);
Apower{1} = eye(n);
for i=2:s+1
    Apower{i} = A*Apower{i-1};
end

B = cell(r+1);
for k=0:r-1
    temp = c(s*k+1)*eye(n);
    for j=1:s-1
        temp = temp + c(s*k+j+1)*Apower{j+1};
    end
    B{k+1} = temp;
end
B{r+1} = c(m+1)*Apower{m-s*r+1};
for j=m-1:-1:s*r
    if j == s*r
       B{r+1} = B{r+1} + c(s*r+1)*eye(n);
    else
       B{r+1} = B{r+1} + c(j+1)*Apower{m-s*r-(m-j)+1};
    end
end

As = Apower{s+1};
P = zeros(n);
for k=r:-1:0
    P = P*As + B{k+1};
end
