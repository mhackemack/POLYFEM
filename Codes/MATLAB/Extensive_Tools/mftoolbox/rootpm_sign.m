function X = rootpm_sign(A,p)
%ROOTPM_SIGN  Matrix Pth root via matrix sign function.
%   X = ROOTPM_SIGN(A,P) computes the principal Pth root X
%   of the matrix A using the matrix sign function and a block
%   companion matrix approach.

if p == 1, X = A; return, end

n = length(A);
podd = rem(p,2);

Y = A;

if podd
    p = 2*p;  % Compensate by squaring at end.
else
    while mod(p,4) == 0
        Y = sqrtm(Y);
        p = round(p/2);
    end
end

if p == 2
    X = sqrtm(Y);
    return
end

% Form C, the block companion matrix.
C = zeros(p*n);
C(end-n+1:end, 1:n) = Y;
C = C + diag(ones(n*(p-1),1),n);

S = signm(C);

X = S(n+1:2*n,1:n);

% Scale factor.
c = 0;
for l = 1:floor(p/4)
    c = c+cos(2*pi*l/p);
end
c = c*4+2;
c = c/p;
X = X/c;
if podd
    X = X*X;
end
