%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          SI-hybrid + MIP matrix functor
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
function P = func_build_SI_hybrid_MIP(lam, input)
% Retrieve Matrices
[L, S] = func_mat_SI_hybrid(lam, input);
[A, B] = func_mat_MIP(lam, input);
% Get Additional Preliminaries
m_quad = input.Quadrature;
m2d = m_quad.moment_to_discrete;
d2m = m_quad.discrete_to_moment;
ndof = input.dof.TotalDoFs;
na = m_quad.NumberAngularDirections;
ntot = na*ndof;
I = eye(ntot);
M = kron(m2d',eye(ndof));
D = kron(d2m,eye(ndof));
% Build Full Matrix
T = L\S; AB = A\B;
MABD = M*AB*D;
P = T + MABD*(T-I);