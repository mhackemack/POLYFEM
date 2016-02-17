%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Iterative Matrix Inversion
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = iterative_matrix_inverse(A,tol)
if nargin == 1, tol = 1e-14; end
I = speye(size(A,1));
% Calculate initial guess
b = eigs(A*A',1);
X0 = (2/b)*A'; norm0 = norm(X0);
% Iterate till convergence
converged = false; iters = 0;
while ~converged
    iters = iters + 1;
    X = X0*(2*I - A*X0);
    % Check convergence
    err = norm(X-X0)/norm0;
    if err < tol, break; end
    X0 = X;
end
% Set outputs
varargout{1} = X;
if nargout > 1, varargout{2} = iters; end