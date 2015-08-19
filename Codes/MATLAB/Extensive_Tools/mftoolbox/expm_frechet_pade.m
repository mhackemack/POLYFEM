function L = expm_frechet_pade(A,E,k)
%EXPM_FRECHET_PADE Frechet derivative of matrix exponential via Pade approx.
%   L = EXPM_FRECHET_PADE(A,E) evaluates the Frechet derivative of
%   the matrix exponential at A in the direction E via scaling and
%   squaring and a Pade approximant of the function tanh(x)/x.
%   L = EXPM_FRECHET_PADE(A,E,k) uses either matrix exponentials
%   (k = 0, the default) or repeated squaring (k = 1) in the final
%   phase of the algorithm.

if nargin < 3, k = 0; end

real_data = isreal(A) && isreal(E);
% Form complex Schur form if A not already upper triangular.
use_Schur = false;
if ~isequal(A,triu(A))
   use_Schur = true;
   [Q,T] = schur(A,'complex'); A = T; E = Q'*E*Q;
end

Abound = 1;
if norm(A,1) <= Abound
   s = 0;
else
   s = ceil( log2(norm(A,1)/Abound) );
end

As = A/2^s;

I = eye(size(A));

m = 8;
% Positive zeros of p8 and q8 in r8 = p8/q8 Pade approximant.
load tau_r8_zeros
% Zeros come in \pm pairs.
a = complex(0, [pzero; -pzero]);
b = complex(0, [qzero; -qzero]);

G = 2^(-s)*E;
for i=1:m
    rhs = (I + As/a(i)) * G + G * (I - As/a(i));
    AA = I + As/b(i); BB = I - As/b(i);
    G = sylvsol(AA, BB, rhs);
end

X = expm(As);
L = (G*X + X*G)/2;
for i=s:-1:1
    if i < s
        if k == 0
           X = expm(2^(-i)*A);
        else
           X = X^2;
        end
    end
    L = X*L + L*X;
end

if use_Schur, L = Q*L*Q'; end
if real_data, L = real(L); end
