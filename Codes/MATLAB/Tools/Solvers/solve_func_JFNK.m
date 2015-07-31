%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve Functor - Custom JFNK
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = solve_func_JFNK(func_list, input, x)
% Strip out necessary function handles
% ------------------------------------
f_rhs = func_list.rhs;
f_resid = func_list.residual;
f_prec = func_list.prec;
% Gather some solver information
% ------------------------------
maxIter = data.solver.maxIterations;
rT = data.solver.relativeTolerance;
aT = data.solver.absoluteTolerance;
F0_norm = norm(f_resid(x,input));
% Loop through maximum number of iterations
for i=1:maxIter
    
    F = f_resid(x,input);
    F_norm = sqrt(F'*F);
    % Exit if converged
    if F_norm <= rT*F_norm + aT, break; end
    % Advance the Newton Step
    
end