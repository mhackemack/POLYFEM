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
function [x, iters] = solve_func_PCG_Eisenstat_Rev1(L, b, x, tol, maxit)
% Check if rhs vector is all zeros
n = length(b);
if ~any(b)
    x = zeros(n,1); iters = 0;
    return
end
% Retrieve necessary matrices
D = (diag(L)); 
b_Ax0 = b - L*x;
DL = tril(L,0);
% Calculate initial values
r = DL\b_Ax0;
rp = D.*r; p = rp;
r0_rp0 = dot(r,rp);
% Loop through iterations
for iters=1:maxit
    % Compute Ap
    t = (DL')\p;
    Ap = t + DL\(p-D.*t);
    % Compute inner products
    r_rp = dot(r,rp);
    a = r_rp/dot(p,Ap);
    % Update solution and residual vectors
    x = x + a*((DL)'\p);
    r = r - a*Ap;
    % Compute new residual inner product norms
    rp = D.*r;
    num = dot(r,rp);
    bb = num/r_rp;
    % Check convergence tolerance
    if abs(num/r0_rp0) < tol
        break;
    end
     p = rp + bb*p;
end