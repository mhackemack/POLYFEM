function [X,Y] = rootpm_schur_newton(A,p)
%ROOTPM_SCHUR_NEWTON  Matrix pth root by Schur-Newton method.
%   [X,Y] = ROOTPM_SCHUR_NEWTON(A,p) computes the principal pth root
%   X of the real matrix A using the Schur-Newton algorithm.
%   It also returns Y = INV(X).

if norm(imag(A),1), error('A must be real.'), end
A = real(A);   % Discard any zero imaginary part.
[Q,R] = schur(A,'real'); % Quasitriangular R.

e = eig(R);
if any (e(find(e == real(e))) < 0 )
   error('A has a negative real eigenvalue: principal pth root not defined')
end

f = factor(p);
k0 = length(find(f == 2)); % Number of factors 2.
q = p/2^(k0);
k1 = k0;
if q > 1

   emax = max(abs(e)); emin = min(abs(e));
   if emax > emin % Avoid log(0).
      k1 = max(k1, ceil( log2( log2(emax/emin) ) ));
   end

   max_arg = norm(angle(e),inf);
   if max_arg > pi/8
      k3 = 1;
      if max_arg > pi/2
         k3 = 3;
      elseif max_arg > pi/4
         k3 = 2;
      end
      k1 = max(k1,k3);
   end

end

for i = 1:k1, R = sqrtm_real(R); end


if q ~= 1
   pw = 2^(-k1); emax = emax^pw; emin = emin^pw;
   if ~any(imag(e))
     % Real eigenvalues.
     if emax > emin
        alpha = emax/emin;
        c = ( (alpha^(1/q)*emax-emin)/( (alpha^(1/q)-1)*(q+1) ) )^(1/q);
     else
        c = emin^(1/q);
     end
   else
     % Complex eigenvalues.
     c = (( emax+emin)/2 )^(1/q);
   end

   X = rootpm_newton(R,q,c);
   for i = 1:k1-k0
       X = X*X;
   end
else
   X = R;
end
Y = Q*(X\Q');  % Return inverse pth root, too.
X = Q*X*Q';
