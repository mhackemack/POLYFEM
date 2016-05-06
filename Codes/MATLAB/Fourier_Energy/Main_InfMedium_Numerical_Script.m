%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Infinite Medium Numerical Simulation Script
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
clear; clc; close all; format long e
% Clear Project Space
% ------------------------------------------------------------------------------
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Home');
inp = 'InfMedium_AllComponents';
addpath(['inputs/',inp]); % This one must be last to properly switch input files
% ------------------------------------------------------------------------------
data = load_user_input();
ng = data.Energy.NumberEngeryGroups;
tol = 1e-8;
maxiters = 1e5;
% ------------------------------------------------------------------------------
nm = data.NumberMaterialsToAnalyze;
JNDD_err = cell(nm, 1); JNDD_NSR = zeros(nm, 1);
JND_err  = cell(nm, 1); JND_NSR  = zeros(nm, 1);
% Loop through materials to analyze and perform analysis
for m=1:nm
    fprintf('Calculating material %d of %d.\n',m,nm)
    % Generate unit source and allocate memory
    Q = rand(ng,1);
    % Get energy bounds from first component in material
    tcompnames = data.Materials{m}.ComponentNames;
    tcompdens = data.Materials{m}.ComponentDensities;
    load([glob.XS_path,tcompnames{1},'/Energy_Bounds.mat']); E = mat;
    Eave = (E(1:ng) + E(2:(ng+1)))./2; Ediff = E(1:ng) - E(2:(ng+1));
    % Zero out total and scattering cross sections
    T = zeros(data.Energy.NumberEngeryGroups);
    S0 = zeros(data.Energy.NumberEngeryGroups);
    % Loop through components within the material
    for c=1:data.Materials{m}.NumberComponents
        % Add total xs contribution
        load([glob.XS_path,tcompnames{c},'/MT_1.mat']);
        T = T + tcompdens(c)*diag(mat);
        % Add P0 scattering xs contribution
        load([glob.XS_path,tcompnames{c},'/MT_2500.mat']);
        S0 = S0 + tcompdens(c)*mat(:,:,1);
    end
    % Solve for exact reference solution
    phiref = (T-S0)\Q;
    % Solve for fast groups
    fg = data.Energy.FastGroups;
    tg = data.Energy.ThermalGroups;
    ntg = length(tg);
    fphi = (T(fg,fg)-S0(fg,fg))\Q(fg,1);
    fsrc = tril(S0(tg,fg)*fphi);
    % Jacobi + NO WGC (DSA + 1G DSA)
    % --------------------------------------------------------------------------
    phi = zeros(ng,1);
    phi(fg,1) = fphi; phi0 = phi;
    err = [];
    % Calculate energy-collapsed XS values
%     F = T(tg,tg)\S0(tg,tg);
%     F = (T(tg,tg) - diag(diag(S0(tg,tg))))\(tril(S0(tg,tg),-1) + triu(S0(tg,tg), 1));
    F = T(tg,tg)\(tril(S0(tg,tg),-1) + triu(S0(tg,tg), 1));
    [V,D] = eig(F); D=(diag(D));
    [~,ind] = max(abs(D));
    V = V(:,ind) / sum(V(:,ind));
    ave_siga = sum((T(tg,tg)-S0(tg,tg))*V);
    % Calculate diffusion xs values
    D = 1./(3*diag(T));
    siga = diag(T) - diag(S0);
    % Loop through iterations
    for i=1:maxiters
        % Sweep
        phi(tg) = T(tg,tg)\(fsrc + Q(tg) + S0(tg,tg)*phi(tg));
        tphi = phi;
        % Loop through thermal groups - perform group-by-group DSA step
        for g=tg(1):tg(end)
            phi(g) = phi(g) + siga(g)\(S0(g,g)*(phi(g) - phi0(g)));
        end
        % Perform EC DSA step
        dp = phi(tg) - phi0(tg);
        rhs = sum(tril(S0(tg,tg),-1)*(dp)) + sum(triu(S0(tg,tg), 1)*(dp));
        dphi = ave_siga\rhs;
        phi(tg) = phi(tg) + V*dphi;
        % Perform convergence testing and update fluxes
        err = [err;max(abs(phi-phi0)./abs(phi))];
        if err(end) < tol
            break;
        end
        phi0 = phi;
    end
    JNDD_err{m} = err;
    JNDD_NSR(m) = JNDD_err{m}(end)/JNDD_err{m}(end-1);
    % Jacobi + NO WGC (1G DSA)
    % --------------------------------------------------------------------------
    phi = zeros(ng,1);
    phi(fg,1) = fphi; phi0 = phi;
    err = [];
    % Calculate energy-collapsed XS values
    F = T(tg,tg)\S0(tg,tg);
    [V,D] = eig(F); D=(diag(D));
    [~,ind] = max(abs(D));
    V = V(:,ind) / sum(V(:,ind));
    ave_siga = sum((T(tg,tg)-S0(tg,tg))*V);
    % Calculate diffusion xs values
    D = 1./(3*diag(T));
    siga = diag(T) - diag(S0);
    % Loop through iterations
    for i=1:maxiters
        % Sweep
        phi(tg) = T(tg,tg)\(fsrc + Q(tg) + S0(tg,tg)*phi(tg));
        tphi = phi;
        % Perform EC DSA step
        rhs = sum(S0(tg,tg)*(phi(tg) - phi0(tg)));
        dphi = ave_siga\rhs;
        phi(tg) = phi(tg) + V*dphi;
        % Perform convergence testing and update fluxes
        err = [err;max(abs(phi-phi0)./abs(phi))];
        if err(end) < tol
            break;
        end
        phi0 = phi;
    end
    JND_err{m} = err;
    JND_NSR(m) = JND_err{m}(end)/JND_err{m}(end-1);
    % GS + WGC (TG)
    % --------------------------------------------------------------------------
%     phi = zeros(ng,1);
%     phi(fg,1) = fphi; phi0 = phi;
%     err = [];
%     % Calculate energy-collapsed XS values
%     F = (T(tg,tg) - tril(S0(tg,tg),0))\triu(S0(tg,tg),1);
%     [V,D] = eig(F); D=(diag(D));
%     [~,ind] = max(abs(D));
%     V = V(:,ind) / sum(V(:,ind));
%     ave_siga = sum((T(tg,tg)-S0(tg,tg))*V);
%     % Calculate diffusion xs values
%     D = 1./(3*diag(T));
%     siga = diag(T) - diag(S0);
%     % Loop through iterations
%     for i=1:maxiters
%         
%     end
end