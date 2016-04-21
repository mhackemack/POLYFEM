%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Multigroup Acceleration Eigenspace Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Project Space
% ------------------------------------------------------------------------------
clc; close all; format short e; clear;
% Load XS Files
% ------------------------------------------------------------------------------
% Thermal groups
ng = 99;
tg = 43:99;
% Energy bounds
load('Energy_Bounds.mat');
E = mat;
Eave  = (E(1:ng) + E(2:end))./2; Eave = Eave(tg);
Ediff = E(1:ng) - E(2:end);     Ediff = Ediff(tg);
clear mat;
% Total xs
load('MT_1.mat');
T = diag(mat(tg));
clear mat;
% Scattering xs
load('MT_2500.mat');
S = mat(tg,tg,1);
clear mat;
% Split scattering xs
SL = tril(S,-1);
SD = diag(diag(S));
SU = triu(S,1);
% Perform Eigenspace Calculations
% ------------------------------------------------------------------------------
hold on;
% Gauss-Seidel + converge inners (TG)
fprintf('-> Gauss-Seidel + converge inners (TG).\n');
A = (T-SL-SD)\(SU);
[V,D] = eig(A); D = diag(D);
[Dmax,ind] = max(abs(D));
Vmax = V(:,ind); Vmax = Vmax / sum(Vmax);
plot(Eave,Vmax./Ediff);
fprintf('  -> Determinant:      %10.5e\n',det(V));
fprintf('  -> Condition #:      %10.5e\n',condest(V));
fprintf('  -> # unique Evalues: %d\n',length(uniquetol(abs(D))));
fprintf('  -> # max Evalue:     %d\n\n',Dmax);
% Gauss-Seidel + do NOT converge inners (MTG)
fprintf('-> Gauss-Seidel + do NOT converge inners (MTG).\n');
A = (T-SL)\(SD+SU);
[V,D] = eig(A); D = diag(D);
[Dmax,ind] = max(abs(D));
Vmax = V(:,ind); Vmax = Vmax / sum(Vmax);
plot(Eave,Vmax./Ediff);
fprintf('  -> Determinant:    %10.5e\n',det(V));
fprintf('  -> Condition #:    %10.5e\n',condest(V));
fprintf('  -> # unique Evalues: %d\n',length(uniquetol(abs(D))));
fprintf('  -> # max Evalue:     %d\n\n',Dmax);
% Jacobi + converge inners (no name yet)
fprintf('-> Jacobi + converge inners (no name yet).\n');
A = (T-SD)\(SL+SU);
[V,D] = eig(A); D = diag(D);
[Dmax,ind] = max(abs(D));
Vmax = V(:,ind); Vmax = Vmax / sum(Vmax);
plot(Eave,Vmax./Ediff);
fprintf('  -> Determinant:    %10.5e\n',det(V));
fprintf('  -> Condition #:    %10.5e\n',condest(V));
fprintf('  -> # unique Evalues: %d\n',length(uniquetol(abs(D))));
fprintf('  -> # max Evalue:     %d\n\n',Dmax);
% Jacobi + do NOT converge inners (MJA - currently used in PDT)
fprintf('-> Jacobi + do NOT converge inners (MJA - currently used in PDT).\n');
A = (T)\(SL+SD+SU);
[V,D] = eig(A); D = diag(D);
[Dmax,ind] = max(abs(D));
Vmax = V(:,ind); Vmax = Vmax / sum(Vmax);
plot(Eave,Vmax./Ediff);
fprintf('  -> Determinant:    %10.5e\n',det(V));
fprintf('  -> Condition #:    %10.5e\n',condest(V));
fprintf('  -> # unique Evalues: %d\n',length(uniquetol(abs(D))));
fprintf('  -> # max Evalue:     %d\n\n',Dmax);
% Complete plotting procedures
box on;
set(gca,'xscale','log');
set(gca,'XGrid','on','XMinorGrid','off');
set(gca,'YGrid','on','YMinorGrid','off');
legend('GS + WGC','GS + NO WGC','Jacobi + WGC','Jacobi + NO WGC');