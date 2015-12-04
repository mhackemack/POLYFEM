%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate Infinite Medium P0 Eigenspectrum
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
function [eig_Jac, eig_GS] = calculate_inf_medium_P0_eigenspectrum(E, T, S)
if size(T,2) == 1, T = diag(T); end
SD0 = tril(squeeze(S(:,:,1))); SU0 = triu(squeeze(S(:,:,1)),1);
% Jacobi iteration matrices
A = T - SD0; B = SU0;
% Two-grid iteration matrices
A = T; B = SD0 + SU0;
