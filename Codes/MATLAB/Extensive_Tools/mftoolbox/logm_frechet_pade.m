function L = logm_frechet_pade(A,E)
%LOGM_FRECHET_PADE Frechet derivative of matrix logarithm via Pade approx.
%   L = LOGM_FRECHET_PADE(A,E) evaluates the Frechet
%   derivative of the matrix logarithm at A in the direction E via the
%   inverse scaling and squaring and a Pade approximant of the function
%   tanh(x)/x.  A must have no eigenvalues on the negative real axis.

real_data = isreal(A) && isreal(E);
% Form complex Schur form if A not already upper triangular.
use_Schur = false;
if ~isequal(A,triu(A))
   use_Schur = true;
   [Q,T] = schur(A,'complex'); A = T; E = Q'*E*Q;
end

I = eye(size(A));
B = A;
for i = 0:inf
    if norm(B-I,1) <= 1-1/exp(1), s = i; break, end
    B = sqrtm(B);
    Aroot{i+1} = B;
end

% Positive zeros of p8 and q8 in r8 = p8/q8 Pade approximant.
load tau_r8_zeros
% Zeros come in \pm pairs.
a = complex(0, [pzero; -pzero]);
b = complex(0, [qzero; -qzero]);

E = 2^s*E;
for i = 1:s
    E = sylvsol(Aroot{i},Aroot{i},E);
end

G = sylvsol(B,B,E);

X = logm(B);

for i=8:-1:1
    rhs = (I + X/b(i)) * G + G * (I - X/b(i));
    AA = I + X/a(i); BB = I - X/a(i);
    G = sylvsol(AA, BB, rhs);
end

L = 2*G;
if use_Schur, L = Q*L*Q'; end
if real_data, L = real(L); end
