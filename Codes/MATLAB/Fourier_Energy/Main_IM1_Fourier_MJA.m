%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          XS Fourier Analysis Script
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
clear; clc; close all; format long e
% Clear Project Space
% ------------------------------------------------------------------------------
fpath = get_path(); addpath(fpath);
global glob; glob = get_globals('Home');
inp = 'IM1_AllComponents';
addpath(['inputs/',inp]); % This one must be last to properly switch input files
% ------------------------------------------------------------------------------
data = load_user_input();
% ------------------------------------------------------------------------------
ny = 4e2; x = linspace(0,4*pi,ny)';
nm = data.NumberMaterialsToAnalyze;
ng = data.Energy.NumberEngeryGroups;
nmom = data.Energy.PnOrder+1;
tg = data.Energy.ThermalGroups;
y_P0_noaccel = zeros(ny,nm);
y_P0_accel   = zeros(ny,nm);
y_P1_noaccel = zeros(ny,nm);
y_P1_accel   = zeros(ny,nm);
% ------------------------------------------------------------------------------
func_P0_noaccel = get_MJA_fourier_func('unaccelerated',0);
func_P0_accel   = get_MJA_fourier_func('accelerated',0);
% func_P1_noaccel = get_MJA_fourier_func('unaccelerated',1);
% func_P1_accel   = get_MJA_fourier_func('accelerated',1);
% ------------------------------------------------------------------------------
mat_str_out = cell(1,data.NumberMaterialsToAnalyze);
P0_noaccel_out = cell(ny+1,data.NumberMaterialsToAnalyze+1);
P0_accel_out = cell(ny+1,data.NumberMaterialsToAnalyze+1);
% P1_noaccel_out = cell(ny+1,data.NumberMaterialsToAnalyze+1);
% P1_accel_out = cell(ny+1,data.NumberMaterialsToAnalyze+1);
P0_noaccel_out{1,1} = 'lambda'; P0_accel_out{1,1} = 'lambda';
% P1_noaccel_out{1,1} = 'lambda'; P1_accel_out{1,1} = 'lambda';
% ------------------------------------------------------------------------------
% Loop through materials to analyze and perform analysis
for m=1:data.NumberMaterialsToAnalyze
    tmat = data.Materials{m};
    matname = data.Materials{m}.MaterialName;
    mat_str_out{m} = matname;
    P0_noaccel_out{1,m+1} = matname; P0_accel_out{1,m+1} = matname;
%     P1_noaccel_out{1,m+1} = matname; P1_accel_out{1,m+1} = matname;
    tcompnames = data.Materials{m}.ComponentNames;
    tcompdens = data.Materials{m}.ComponentDensities;
    % Zero out total and scattering cross sections
    T = zeros(ng,ng); S = zeros(ng,ng,nmom);
    % Loop through components within the material
    for c=1:data.Materials{m}.NumberComponents
        % Add total xs contribution
        load([glob.XS_path,tcompnames{1},'/MT_1.mat']);
        T = T + tcompdens(c)*diag(mat);
        % Add P0 scattering xs contribution
        load([glob.XS_path,tcompnames{1},'/MT_2500.mat']);
        S = S + tcompdens(c)*mat;
    end
    % Form TG spectrum and diffusion coefficients
    A = T(tg,tg)\S(tg,tg,1);
    [V,D] = eig(A); D = diag(D);
    [eval,Ei] = max(abs(D)); D = [];
    V = V(:,Ei); V = V / sum(V); Vtg = V;
    D0 = (1/3)./diag(T); D1 = zeros(ng,1);
    for g=1:ng
        D1(g) = 1/(3*(T(g,g) - sum(S(:,g,2))));
    end
    % Generate fourier mode distributions
    display(sprintf('Computing Fourier modes for material %d of %d.',m,data.NumberMaterialsToAnalyze))
    for i=1:ny
%         display(sprintf('    Computing mode: %d of %d',i,ny))
        y_P0_noaccel(i,m) = func_P0_noaccel(x(i), T(tg,tg,1), S(tg,tg,1));
        y_P0_accel(i,m)   = func_P0_accel(x(i), T(tg,tg,1), S(tg,tg,1), D0(tg), V);
%         y_P1_noaccel(i,m) = func_P1_noaccel(x(i), T(tg,tg,1), S(tg,tg,1));
%         y_P1_accel(i,m)   = func_P1_accel(x(i), T(tg,tg,1), S(tg,tg,1), D1(tg), V);
        % Place fourier results into global data structures
        P0_noaccel_out{1+i,m+1} = num2str(y_P0_noaccel(i,m),'%16.8e');
        P0_accel_out{1+i,m+1}   = num2str(y_P0_accel(i,m),'%16.8e');
%         P1_noaccel_out{1+i,m+1} = num2str(y_P1_noaccel(i,m),'%16.8e');
%         P1_accel_out{1+i,m+1}   = num2str(y_P1_accel(i,m),'%16.8e');
    end
end
% ------------------------------------------------------------------------------
for i=1:ny
    P0_noaccel_out{1+i,1} = num2str(x(i),'%16.8e');
    P0_accel_out{1+i,1}   = num2str(x(i),'%16.8e');
%     P1_noaccel_out{1+i,1} = num2str(x(i),'%16.8e');
%     P1_accel_out{1+i,1}   = num2str(x(i),'%16.8e');
end
% ------------------------------------------------------------------------------
% Check if output directory exists
if ~isequal(exist([glob.output_path,data.OutputName], 'dir'),7),mkdir([glob.output_path,data.OutputName]); end
% Print output files
cell2csv([glob.output_path,data.OutputName,'/Fourier_MJA_P0_noaccel'],P0_noaccel_out);
cell2csv([glob.output_path,data.OutputName,'/Fourier_MJA_P0_accel'],P0_accel_out);
% cell2csv([glob.output_path,data.OutputName,'/Fourier_MJA_P1_noaccel'],P1_noaccel_out);
% cell2csv([glob.output_path,data.OutputName,'/Fourier_MJA_P1_accel'],P1_accel_out);
