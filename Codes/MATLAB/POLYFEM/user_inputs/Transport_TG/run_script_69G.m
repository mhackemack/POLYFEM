%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          69G Graphite Two-Grid Run Script
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
% Clear Project Space
% ------------------------------------------------------------------------------
clear; clc; close all; format long e;
% Define some group/order information
% ------------------------------------------------------------------------------
% dir = '119G_graphite';
% ng = 119; Pn = 0;
% fg = 1:62; nfg = length(fg);
% tg = 63:119; ntg = length(tg);
dir = '69G_graphite';
ng = 69; Pn = 0;
fg = 1:28; nfg = length(fg);
tg = 29:69; ntg = length(tg);
% Retrieve XS Data
% ------------------------------------------------------------------------------
% Get Energy Bounds
load([dir,'/Energy_Bounds.mat']);
E = mat; clear mat;
Ediff = E(1:end-1) - E(2:end);
Eave = (E(1:end-1) + E(2:end))./2;
% Get Total Cross Sections
load([dir,'/MT_1.mat']);
txs = mat; clear mat;
T = diag(txs); Tave = txs./Ediff;
% Get Transfer Cross Sections
load([dir,'/MT_2500.mat']);
S = mat; clear mat;
S = S(:,:,1:Pn+1);
% Calculate Infinite Medium Eigenvalue
% ------------------------------------------------------------------------------
A = (T - tril(S(:,:,1)))\triu(S(:,:,1),1);
% A = (diag(txs(tg)) - tril(sxs(tg,tg,1)))\triu(sxs(tg,tg,1),1);
[V,D] = eig(A); D = diag(D);
[eval,Ei] = max(abs(D));
V = V(:,Ei); V = V / sum(V); Vtg = V(tg);
if Pn==0
    D = (1/3)./txs;
    Dave = sum(D.*V);
    siga = 0;
    for g=1:ng
        siga = siga + V(g)*txs(g);
        for gg=1:ng
            siga = siga - S(g,gg,1)*V(gg);
        end
    end
elseif Pn==1
    D = zeros(ng,1); siga = 0;
    for g=1:ng
        siga = siga + V(g)*txs(g); tt = 0;
        for gg=1:ng
            tt = tt + S(gg,g,2);
            siga = siga - S(g,gg,1)*V(gg);
        end
        D(g) = 1/(3*(txs(g) - tt));
    end
    Dave = sum(D.*V);
end
% P0 Fourier Analysis
% ------------------------------------------------------------------------------
noaccel_func_P0 = get_2G_fourier_func('unaccelerated',0);
accel_func_P0 = get_2G_fourier_func('accelerated',0);
n = 1e2; x = linspace(1e-8,2,n);
y_P0_noaccel = zeros(n,1); y_P0_accel = zeros(n,1);
for i=1:n
    fprintf('P0 Fourier iterate: %d of %d\n',i,n);
    y_P0_noaccel(i) = noaccel_func_P0(x(i), T, S);
    y_P0_accel(i) = accel_func_P0(x(i), T, S, D, V);
end
% P1 Fourier Analysis
% ------------------------------------------------------------------------------
% if Pn == 1
%     noaccel_func_P1 = get_2G_fourier_func('unaccelerated',1);
%     accel_func_P1 = get_2G_fourier_func('accelerated',1);
%     n = 2e2; x = linspace(1e-6,2,n)';
%     y_P1_noaccel = zeros(n, 1); y_P1_accel = zeros(n, 1);
%     for i=1:n
%         fprintf('P1 Fourier iterate: %d of %d\n',i,n);
%         y_P1_noaccel(i) = noaccel_func_P1(x(i), T, S);
%     end
% end
% P0 Analytical Analysis
% ------------------------------------------------------------------------------
A = T - S(:,:,1); b = [1;zeros(ng-1,1)];
sol_ana = A\b;
% P0 Numerical Analysis - no acceleration
% ------------------------------------------------------------------------------
q = [1;zeros(ng-1,1)]; sol = zeros(ng,1);
itmax = 1e8; tol = 1e-8;
% Loop through the fast groups
for g=fg(1):fg(end)
    b = q(g); A = txs(g) - S(g,g,1);
    b = b + S(g,:,1)*sol;
    sol(g) = A\b;
end
sol0 = sol; oldgnorm = 1;
% Perform iterations over thermal groups
SR_noaccel = [];
for m=1:itmax
    % Loop through thermal groups
    for g=tg(1):tg(end)
        b = q(g); A = txs(g);
        b = b + S(g,:,1)*sol;
        sol(g) = A\b;
    end
    % Calculate error
    inf_norm = max(abs(sol-sol0));
    gnorm = sum((sol-sol0).^2);
    SR_noaccel = [SR_noaccel; sqrt(gnorm/oldgnorm)];
    err = gnorm / sum(sol.^2);
    fprintf('Iteration # %d: %0.6e\n', m, inf_norm);
    if inf_norm < tol, break; end
    % Update iterative values
    sol0 = sol; oldgnorm = gnorm;
end
sol_num_noaccel = sol;
% P0 Numerical Analysis - with acceleration
% ------------------------------------------------------------------------------
q = [1;zeros(ng-1,1)]; sol = zeros(ng,1);
% Loop through the fast groups
for g=fg(1):fg(end)
    b = q(g); A = txs(g) - S(g,g,1);
    b = b + S(g,:,1)*sol;
    sol(g) = A\b;
end
sol0 = sol; oldgnorm = 1;
% Perform iterations over thermal groups
SR_accel = [];
for m=1:itmax
    % Loop through thermal groups
    for g=tg(1):tg(end)
        b = q(g); A = txs(g);
        b = b + S(g,:,1)*sol;
        sol(g) = A\b;
    end
    % Perform acceleration
    R = 0;
    for g=tg(1):tg(end)
        for gg=g+1:tg(end)
            R = R + S(g,gg)*(sol(gg)-sol0(gg));
        end
    end
    dx = R / siga;
    sol = sol + dx*V;
    % Calculate error
    inf_norm = max(abs(sol-sol0));
    gnorm = sum((sol-sol0).^2);
    SR_accel = [SR_accel; sqrt(gnorm/oldgnorm)];
    err = gnorm / sum(sol.^2);
    fprintf('Iteration # %d: %0.6e\n', m, inf_norm);
    if inf_norm < tol, break; end
    % Update iterative values
    sol0 = sol; oldgnorm = gnorm;
end
sol_num_accel = sol;
