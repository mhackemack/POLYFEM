function X = rootpm_real(A,p)
%ROOTPM_REAL  Pth root of real matrix via real Schur form.
%   X = ROOTPM_REAL(A,P) is a Pth root of the nonsingular matrix A
%   (A = X^P).  When A has no eigenvalues on the negative real axis
%   then X is the principal pth root, which is real when A is real.

%   Implementation is complicated by use of zero (and even "-1")
%   subscripts in algorithm.  Dealt with by increasing "k" by 1
%   each time and moving "-1" cases in front of loops.

if norm(imag(A),1), error('A must be real.'), end
if p == 1, X = A; return; end

n = length(A);
[Q,R] = schur(A,'real');

% m blocks: i'th has order s(i), starts at t(i).
[m,s,t] = quasitriang_struct(R);
U = zeros(n); B = cell(p); V = cell(p-1);
V{1} = eye(n);  % U^k == V{k+1}.

for j=1:m
    rj = t(j):t(j)+s(j)-1;
    if s(j) == 1
       % If A is real, X will be real unless R(p,p)<0 on the next line.
       U(rj,rj) = R(rj,rj)^(1/p);
    else
       U(rj,rj) = root_block(R(rj,rj),p);
    end
    for k = 0:p-2, kk = k+1; V{kk}(rj,rj) = U(rj,rj)^(k+1); end
    for i=j-1:-1:1
        ri = t(i):t(i)+s(i)-1;
        for k = 0:p-2
            kk = k+1;
            B{kk} = zeros(s(i),s(j));
            for ell = i+1:j-1
                rell = t(ell):t(ell)+s(ell)-1;
                B{kk} = B{kk} + U(ri,rell)*V{kk}(rell,rj);
            end
        end
        rhs = R(ri,rj) - B{p-2 +1};
        for k = 0:p-3
            rhs = rhs - V{p-3-k +1}(ri,ri)*B{k+1};
        end
        coeff = kron( eye(s(j)), V{p-2 +1}(ri,ri) ) ...
               + kron( V{p-2 +1}(rj,rj).', eye(s(i)) );
        for k = 1:p-2
            coeff = coeff + kron( V{k-1 +1}(rj,rj).', V{p-2-k +1}(ri,ri) );
        end
        y = coeff\rhs(:);
        rhs(:) = y;                  % `Un-vec' the solution.
        U(ri,rj) = rhs;
        for k = 0:p-2
            if k == 0
               S = U(ri,rj);
            else
               S = V{k-1 +1}(ri,ri)*U(ri,rj) + U(ri,rj)*V{k-1 +1}(rj,rj) ...
                   + B{k-1 +1};
               for ell = 1:k-1
                   S = S + V{k-ell-1 +1}(ri,ri)*U(ri,rj)*V{ell-1 +1}(rj,rj);
               end
               for ell = 0:k-2
                   S = S + V{k-2-ell +1}(ri,ri)*B{ell +1};
               end
            end
            V{k+1}(ri,rj) = S;
        end
   end
end

X = Q*U*Q';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = root_block(R,p)
%ROOT_BLOCK  Pth root of a real 2x2 matrix with complex conjugate eigenvalues.

if norm(imag(R)) ~= 0 || any(size(R) - [2,2])
   error('Matrix must be real, of dimension 2.')
end

r11 = R(1,1); r12 = R(1,2);
r21 = R(2,1); r22 = R(2,2);

theta = (r11 + r22) / 2;
musq = (-(r11 - r22)^2 - 4*r21*r12) / 4;
mu = sqrt(musq);

if musq <= 0
   error('Matrix must have non-real complex conjugate eigenvalues.')
end

r = sqrt(theta^2+musq);
phi = angle(complex(theta,mu));
rootp = r^(1/p)*exp(i*phi/p);

alpha = real(rootp);
beta = imag(rootp);

X = alpha*eye(2) + (beta/mu)*(R - theta*eye(2));
