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
function x = solve_func_PCG_Eisenstat(L, b, x, tol, maxit)
% Check if rhs vector is all zeros
n = length(b);
if ~any(b)
    x = zeros(n,1);
    return
end
% Retrieve necessary matrices
D = (diag(L));
b_Ax0 = b - L*x;
L = tril(L,-1); 
% Calculate initial values
r = (diag(D)+L)\b_Ax0; D = full(D);
rp = D.*r; p = rp;
r0_rp0 = dot(r,rp);
% Loop through iterations
for iter=1:maxit
    Ap = compute_Ap(D, L, p);
    r_rp = dot(r,rp);
    denom = dot(p,Ap);
    a = r_rp/denom;
    
    tmp = (diag(D) + L)'\p;
%     tmp = backward_substitution(D, L, p);
    x = x + a*tmp;
    r = r - a*Ap;
    
    rp = D.*r;
    num = dot(r,rp);
    bb = num/r_rp;
    % Check convergence tolerance
    if abs(num/r0_rp0) < tol
        break;
    end
     p = rp + bb*p;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = residual_norm(x,y)
out = sqrt(sum(x.*y));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = forward_substitution(D, L, b)
n = length(b); out = b;
out(1) = b(1)/D(1,1);
% Loop forwards through vector length
for i=2:n
    c = L(i,:); [~,col] = find(c);
    for j=1:length(col)
        out(i) = out(i) - c(col(j))*out(col(j));
    end
    out(i) = out(i) / D(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = backward_substitution(D, L, b)
n = length(b); out = b;
% Loop backwards through vector length
for i=n:-1:1
    out(i) = b(i);
    c = L(i,:); [~,col] = find(c);
    tcol = n-col+1;
    for j=1:length(col)
        out(i) = out(i) - c(col(j))*out(tcol(j));
    end
    out(i) = out(i) / D(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = compute_Ap(D, L, p)
t = (diag(D) + L)'\p;
% t = backward_substitution(D, L, p);
out = t + forward_substitution(D, L, p - D.*t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%