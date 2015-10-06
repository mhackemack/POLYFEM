
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get appropriate two-grid fourier function handle
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
function out = get_2G_fourier_func(ftype, Pn)

if Pn == 0
    if strcmpi(ftype, 'unaccelerated')
        out = @P0_unaccel_func;
    elseif strcmpi(ftype, 'accelerated')
        out = @P0_accel_func;
    else
        error('Do not know what you mean...');
    end
elseif Pn == 1
    if strcmpi(ftype, 'unaccelerated')
        out = @P1_unaccel_func;
    elseif strcmpi(ftype, 'accelerated')
        out = @P1_accel_func;
    else
        error('Do not know what you mean...');
    end
else
    error('Only P0 and P1 right now.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fourier Function Handles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = P0_unaccel_func(val, T, S)
ng = size(T,1); I = eye(ng);
Sd = tril(S(:,:,1),0); Su = triu(S(:,:,1),1);
if abs(val) < 1e-10
    Tmat = diag(1./diag(T));
else
    Tmat = diag(atan(val./diag(T)))/val;
end
L1 = I - Tmat*Sd;
L2 = Tmat*Su;
out = max(abs(eig(L1\L2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = P0_accel_func(val, T, S, D, E)
ng = size(T,1); I = eye(ng);
Sd = tril(S(:,:,1),0); Su = triu(S(:,:,1),1);
D = diag(D);
if abs(val) < 1e-10
    Tmat = diag(1./diag(T));
else
    Tmat = diag(atan(val./diag(T)))/val;
end
L1 = I - Tmat*Sd; L2 = Tmat*Su;
L3 = E*(sum((val^2*D+T-Sd-Su)*E))^(-1)*sum(Su,1);
out = max(abs(eig((I+L3)*(L1\L2) - L3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = P1_unaccel_func(val, T, S)
ng = size(T,1); I = eye(ng);
Sd0 = tril(S(:,:,1)); Su0 = triu(S(:,:,1),1);
Sd1 = tril(S(:,:,2)); Su1 = triu(S(:,:,2),1);
if abs(val) < 1e-10
    error('Cannot limit my way out of this one...');
else
    Tmat = diag(atan(val./diag(T)))/val;
    TTmat = T*diag(atan(val./diag(T)))/val;
end
L1 = [I-Tmat*Sd0, 3*1i/val*(I-TTmat)*Sd1;...
      1i/val*(I-TTmat)*Sd0, I-3/val^2*T*(I-TTmat)*Sd1];
L2 = [Tmat*Su0, -3*1i/val*(I-TTmat)*Su1;...
      -1i/val*(I-TTmat)*Sd0, 3/val^2*T*(I-TTmat)*Su1];
out = max(abs(eig(L1\L2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = P1_accel_func(val, T, S, D, E)
ng = size(T,1); I = eye(ng); z = zeros(ng);
Sd0 = tril(S(:,:,1)); Su0 = triu(S(:,:,1),1);
Sd1 = tril(S(:,:,2)); Su1 = triu(S(:,:,2),1);
D = diag(D);
if abs(val) < 1e-10
    error('Cannot limit my way out of this one...');
else
    Tmat = diag(atan(val./diag(T)))/val;
    TTmat = T*diag(atan(val./diag(T)))/val;
end
L1 = [I-Tmat*Sd0, 3*1i/val*(I-TTmat)*Sd1;...
      1i/val*(I-TTmat)*Sd0, I-3/val^2*T*(I-TTmat)*Sd1];
L2 = [Tmat*Su0, -3*1i/val*(I-TTmat)*Su1;...
      -1i/val*(I-TTmat)*Sd0, 3/val^2*T*(I-TTmat)*Su1];
L3 = [E*(sum((val^2*D+T-Sd0-Su0)*E))^(-1)*sum(Su0,1), z; z, z];
LL = L1\L2;
out = LL-L3*(eye(2*ng) - LL);
out = max(abs(eig(out)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%