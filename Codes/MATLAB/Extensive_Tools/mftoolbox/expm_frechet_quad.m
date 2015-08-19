function R = expm_frechet_quad(A,E,theta,rule,k)
%EXPM_FRECHET_QUAD Frechet derivative of matrix exponential via quadrature.
%   L = EXPM_FRECHET_QUAD(A,E,THETA,RULE) is an approximation to the
%   Frechet derivative of the matrix exponential at A in the direction E
%   intended to have norm of the correct order of magnitude.
%   It is obtained from the repeated trapezium rule (RULE = 'T'),
%   the repeated Simpson rule (RULE = 'S', default),
%   or the repeated midpoint rule (RULE = 'M').
%   L = EXPM_FRECHET_QUAD(A,E,THETA,RULE,k) uses either matrix
%   exponentials (k = 0, the default) or repeated squaring (k = 1)
%   in the final phase of the algorithm.

%   A is scaled so that norm(A/2^s) <= THETA.  Defalt: THETA = 1/2.

if nargin < 3 || isempty(theta), theta = 1/2; end
if nargin < 4 || isempty(rule), rule = 'S'; end
if nargin < 5, k = 0; end

s = ceil( log2(norm(A,1)/theta) );
As = A/2^s;

X = expm(As);

switch upper(rule)

    case 'T'
       R = 2^(-s) * (X*E + E*X)/2 ;

    case 'S'
       Xmid = expm(As/2);
       R = 2^(-s) * (X*E + 4*Xmid*E*Xmid + E*X)/6;

    case 'M'
       Xmid = expm(As/2);
       R = 2^(-s) * Xmid*E*Xmid;

    otherwise
        error('Illegal value of RULE.')

end

for i = s:-1:1
    if i < s
        if k == 0
           X = expm(2^(-i)*A);
        else
           X = X^2;
        end
    end
    R = X*R + R*X;
end
