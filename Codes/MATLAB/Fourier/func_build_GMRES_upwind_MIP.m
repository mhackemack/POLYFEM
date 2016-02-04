%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          GMRES-upwind + MIP matrix functor
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
function P = func_build_GMRES_upwind_MIP(lam, input)
% Retrieve Matrices
T = func_mat_SI_upwind(lam, input);
[A, B] = func_mat_MIP(lam, input);
% Get Additional Preliminaries
I = eye(size(T,1));
% Build Full Matrix
P = I - (T + (A\B)*(T - I));
