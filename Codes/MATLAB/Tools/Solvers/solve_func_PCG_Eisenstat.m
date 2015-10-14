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
function x = solve_func_PCG_Eisenstat(A, b, x, data, varargin)

% Retrieve Matrix and RHS Vector
% A = Lfunc(x, data, varargin);
% b = Rfunc(x, data, varargin);
% Form Operation Matrices and Compute Initial Vectors
D = diag(diag(A)); DL = tril(A);
I = speye(size(D,1)); DLinv = DL\I; DLinvt = DLinv';
AA = DLinv*A*DLinvt*DL';
r = DLinv*(b-A*x); p = D*r; rp = p;
norm0 = residual_norm(r,rp);
% Run through iterations
rT = data.solver.relativeTolerance;
for i=1:data.solver.maxIterations
    % Get new step length
    a = residual_norm(r,rp) / residual_norm(p,AA*p);
    x = x + a*DLinvt*p;
    % Form new residual
    r = r + a*AA*p;
    rp = D*r;
    b = residual_norm(r,rp) / norm0;
    if b < rT, break; end
    % Update search direction
    p = rp + b*p;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = residual_norm(x,y)
out = sqrt(sum(x.*y));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%