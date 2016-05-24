%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get appropriate MTG fourier function handle
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
function out = get_MTG_fourier_func(ftype, Pn)
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
L1 = I - Tmat*Sl;
L2 = Tmat*(Sd + Su);
out = max(abs(eig(L1\L2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = P0_accel_func(val, T, S, D, E)
ng = size(T,1); I = eye(ng);
Sl = tril(S(:,:,1),-1);
Sd = diag(diag((S(:,:,1))));
Su = triu(S(:,:,1),1);
D = diag(D);
if abs(val) < 1e-10
    Tmat = diag(1./diag(T));
else
    Tmat = diag(atan(val./diag(T)))/val;
end
L1 = I - Tmat*Sl;
L2 = Tmat*(Sd + Su);
L3 = E*(sum((val^2*D+T-Sl-Sd-Su)*E))^(-1)*sum((Sd + Su),1);
out = max(abs(eig((I+L3)*(L1\L2) - L3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%