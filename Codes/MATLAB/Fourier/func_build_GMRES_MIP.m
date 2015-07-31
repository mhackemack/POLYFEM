%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          GMRES + MIP matrix functor
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
function P = func_build_GMRES_MIP(lam, input)
% Retrieve Matrices
T = func_mat_GMRES(lam, input);
[A, B] = func_mat_MIP(lam, input);
% Get Additional Preliminaries
I = eye(input.dof.TotalDoFs);
% Build Full Matrix
TT = I-T;
P = TT + (A\B)*(TT - I);