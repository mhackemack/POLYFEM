function [X,k,m] = logm_iss(A)
%LOGM_ISS Matrix logarithm by inverse scaling and squaring method.
%   X = LOGM_ISS(A) computes the logarithm of A, for a matrix with no
%   nonpositive real eigenvalues, using the inverse scaling and squaring
%   method with Pade approximation.
%   Matrix square roots are computed by the product form of the
%   Denman-Beavers teration.
%   [X,K,M] = LOGM_ISS(A) returns the number K of square roots
%   computed and the degree M of the Pade approximant.

e = eig(A);
if any( imag(e) == 0 & real(e) <= 0 )
   error('A must not have any nonpositive real eigenvalues!')
end

n = length(A);

load log_pade_err_opt  % mmax-by-3 matrix DATA.
% mvals = data(:,1);
xvals = data(:,2);

X = A;
k = 0; p = 0; itk = 5;

while 1

    normdiff = norm(X-eye(n),1);
    if normdiff <= xvals(16)

       p = p+1;
       j1 = find(normdiff <= xvals(3:16));
       j1 = j1(1) + 2;
       j2 = find(normdiff/2 <= xvals(3:16));
       j2 = j2(1) + 2;
       if 2*(j1-j2)/3 < itk || p == 2, m = j1; break, end

    end

    [X,M,itk] = sqrtm_dbp(X,1); k = k+1;

end

X = 2^k*logm_pade_pf(X-eye(n),m);
if isreal(A), X = real(X); end
