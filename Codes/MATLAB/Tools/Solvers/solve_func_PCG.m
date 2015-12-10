%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve Functor - Custom PCG
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, iters] = solve_func_PCG(A, b, x, M1, M2, tol, maxiter)
n = length(b); normb = norm(b);
% Check for all zero rhs vector
if normb == 0
    x = zeros(n,1); iters = 0;
    return
end
% Check if preconditioner matrices exist
P1_bool = false; P2_bool = false;
if ~isempty(M1), P1_bool = true; end
if ~isempty(M2), P2_bool = true; end
% Get initial values and setup system
tolb = tol*norm(b); converged = false;
r = b - A*x; normr = norm(r);
% Check if initial guess is close enough 
if normr <= tolb
    iters = 0;
    return
end
% Loop through iterations
iters = 0; rho_1 = 0;
while iters <= maxiter && ~converged
    iters = iters + 1;
    % Appply first preconditioner
    if P1_bool
        y = M1\r;
    else
        y = r;
    end
    % Apply second preconditioner
    if P2_bool
        z = M2\y;
    else
        z = y;
    end
    % Calculate the new search direction
    rho = r'*z;
    if iters == 1
        p = z;
    else
        p = z + (rho / rho_1)*p;
    end
    % Perform matrix action on the krylov vector
    q = A*p;
    pq = p'*q;
    % Form new step length
    aa = rho / pq;
    % Update solution and residual
    x = x + aa*p;
    r = r - aa*q;
    normr = norm(r);
    % Check for convergence
    if normr <= tolb
        converged = true;
    else
        rho_1 = rho;
    end
end