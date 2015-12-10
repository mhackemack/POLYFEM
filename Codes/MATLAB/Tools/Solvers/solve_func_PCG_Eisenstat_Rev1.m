%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve Functor - Custom PCG with Eisenstat Trick
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = solve_func_PCG_Eisenstat_Rev1(L, b, x, tol, maxit)
% Check if rhs vector is all zeros
n = length(b);
if ~any(b)
    x = zeros(n,1);
    return
end
% Retrieve necessary matrices
D = (diag(L)); 
b_Ax0 = b - L*x;
L = tril(L,-1); D = sparse(1:n,1:n,D);
% Calculate initial values
r = (D+L)\b_Ax0;
rp = D.*r; p = rp;
r0_rp0 = dot(r,rp);
% Loop through iterations
for iter=1:maxit
    
end