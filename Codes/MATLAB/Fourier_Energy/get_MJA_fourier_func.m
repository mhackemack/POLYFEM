%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get appropriate MJA fourier function handle
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
function out = get_MJA_fourier_func(ftype, Pn)
if Pn == 0
    if strcmpi(ftype, 'unaccelerated')
        out = @P0_unaccel_func;
    elseif strcmpi(ftype, 'accelerated')
        out = @P0_accel_func;
    else
        error('Do not know what you mean...');
    end
else
    error('Only P0 right now.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fourier Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = P0_unaccel_func(val, T, S)
ng = size(T,1); I = eye(ng);
Sl = tril(S(:,:,1),-1);
Sd = diag(diag((S(:,:,1))));
Su = triu(S(:,:,1),1);
if abs(val) < 1e-10
    Tmat = diag(1./diag(T));
else
    Tmat = diag(atan(val./diag(T)))/val;
end
L1 = Tmat;
% L1 = I - Tmat*(Sl + Sd + Su);
L2 = Tmat*(Sl + Sd + Su);
out = max(abs(eig(L2)));
% out = max(abs(eig(L1\L2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = P0_accel_func(val, T, S, D, E)
ng = size(T,1); I = eye(ng);
% Sl = tril(S(:,:,1),-1);
% Sd = diag(diag((S(:,:,1))));
% Su = triu(S(:,:,1),1);
D = diag(D);
if abs(val) < 1e-10
    Tmat = diag(1./diag(T));
else
    Tmat = diag(atan(val./diag(T)))/val;
end
L1 = Tmat*S;
L2 = E*(sum((val^2*D+T-S)*E))^(-1)*sum(S*(L1-I),1);
out = max(abs(eig(L1 + L2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%