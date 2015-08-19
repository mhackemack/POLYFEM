function X = sqrtm_real(A)
%SQRTM_REAL Square root of real matrix by real Schur method.
%   X = SQRTM_REAL(A) is the principal square root of the real matrix A,
%   computed in real arithmetic by the real Schur method.
%   X is real unless A has a real negative eigenvalue
%   (in this case a real primary square root does not exist).

if norm(imag(A),1), error('A must be real.'), end
A = real(A);   % Discard any zero imaginary part.
n = length(A);
[Q,R] = schur(A,'real'); % Quasitriangular R.

% m blocks: i'th has order s(i), starts at t(i).
[m, s, k] = quasitriang_struct(R);
T = zeros(n);

for j=1:m
    p = k(j):k(j)+s(j)-1;
    if s(j) == 1
       % If A is real, X will be real unless R(p,p)<0 on the next line.
       T(p,p) = sqrt(R(p,p));
    else
       T(p,p) = rsqrt2(R(p,p));
    end
    for r=j-1:-1:1
        rind = k(r):k(r)+s(r)-1;
        rj = k(r+1):k(j)-1;
        if ~isempty(rj)
           prod = T(rind,rj)*T(rj,p);  % Gives [] when rj = [].
        else
           prod = zeros(s(r),s(j));
        end
        B = R(rind,p) - prod;
        % NB Unconjugated transpose on next line for complex case.
        A = kron( eye(s(j)), T(rind,rind) ) + kron( T(p,p).', eye(s(r)) );
        y = A\B(:);
        B(:) = y;                      % `Un-vec' the solution.
        T(rind,p) = B;
   end
end

X = Q*T*Q';

%%%%%%%%%%%%%%%%%%%%%%
function X = rsqrt2(R)
%RSQRT2  Real square root of a real 2x2 matrix with complex conjugate
%        eigenvalues.

r11 = R(1,1); r12 = R(1,2);
r21 = R(2,1); r22 = R(2,2);

theta = (r11 + r22) / 2;
musq = (-(r11 - r22)^2 - 4*r21*r12) / 4;
mu = sqrt(musq);

if theta > 0
   alpha = sqrt( (theta + sqrt(theta^2+musq))/2 );
else
   alpha = mu / sqrt( 2*(-theta + sqrt(theta^2+musq)) );
end

X = (alpha-theta/(2*alpha)) * eye(2) + R/(2*alpha);
