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
function [PI_eig_Jac,PI_eig_GS,PI_eig_MTG,eig_Jac,eig_GS,eig_MTG] = calculate_inf_medium_P0_eigenspectrum(g, T, S)
if size(T,2) == 1, T = diag(T); end
SD0 = tril(squeeze(S(:,:,1))); SU0 = triu(squeeze(S(:,:,1)),1);
% Gauss-Seidel iteration matrices - MATLAB eig
A = T(g,g) - SD0(g,g); B = SU0(g,g); C = A\B;
[eig_GS.EigenSpectrum,De] = eig(C);
eig_GS.EigenValue = diag(De);
% Gauss-Seidel iteration matrices - Power Iteration
[PI_eig_GS.EigenSpectrum,PI_eig_GS.EigenValue,it] = power_method(C, ones(length(g),1),1e4,1e-12);
% Jacobi iteration matrices - MATLAB eig
A = T(g,g); B = SD0(g,g) + SU0(g,g); C = A\B;
[eig_Jac.EigenSpectrum,De] = eig(C);
eig_Jac.EigenValue = diag(De);
% Jacobi iteration matrices - Power Iteration
[PI_eig_Jac.EigenSpectrum,PI_eig_Jac.EigenValue,it] = power_method(C, ones(length(g),1),1e4,1e-12);
% MTG iteration matrices - MATLAB eig
A = T(g,g) - tril(S(g,g,1),-1); B = diag(diag(S(g,g,1))) + SU0(g,g); C = A\B;
[eig_MTG.EigenSpectrum,De] = eig(C);
eig_MTG.EigenValue = diag(De);
% MTG iteration matrices - Power Iteration
[PI_eig_MTG.EigenSpectrum,PI_eig_MTG.EigenValue,it] = power_method(C, ones(length(g),1),1e4,1e-12);