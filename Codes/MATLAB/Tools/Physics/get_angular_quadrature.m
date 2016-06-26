%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate Angular Quadrature Set
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_angular_quadrature(data, dim)
% Get Angular Quadrature Data
% ---------------------------
AQName = data.QuadType;
% Get Single Octant Data
% ----------------------
if strcmpi(AQName, 'Manual')
    x = data.QuadAngles; nx = size(x,1);
    x = [x,zeros(nx,3-dim)];
    w = data.QuadWeights;
    adims = [1,2,3];
else
    if dim == 1
        [x,w] = lgwt(data.SnLevels,-1,1);
        w = w*2;
        adims = [1,2,3];
    else
        % Form angle directions based on polar angle - defaults to z
        if isfield(data, 'PolarDimension')
            pdim = data.PolarDimension;
            if pdim == 1
                adims = [2,3,1];
            elseif pdim == 2
                adims = [1,3,2];
            elseif pdim == 3
                adims = [1,2,3];
            end
        else
            adims = [1,2,3];
        end
        % Switch based on quadrature type
        if strcmp(AQName, 'LDFE')
            [x,w] = get_LDFE_quad(dim, data.SnLevels, adims);
        elseif strcmp(AQName, 'GLC')
            [x,w] = get_GLC_quad(dim, data.SnLevels, adims);
        elseif strcmp(AQName, 'PGLC')
            [x,w] = get_PGLC_quad(dim, data.PolarLevels, data.AzimuthalLevels, adims);
        elseif strcmp(AQName, 'LS')
            [x,w] = get_LS_quad(dim, data.SnLevels, adims);
        elseif strcmp(AQName, 'TriGLC')
            [x,w] = get_TriGLC_quad(dim, data.PolarLevels, data.AzimuthalLevels, adims);
        else
            error('Cannot determine angular quadrature type.')
        end
    end
end
% Specify Angular Quadrature Norm
if strcmpi(AQName, 'Manual')
    a_norm = sum(w);
else
    if dim == 1
        a_norm = 2;
    else
        a_norm = 2^(dim-1)*pi;
    end
end
% Calculate remaining information
if strcmpi(AQName, 'Manual')
    opp_ind = determine_opposite_angle(x);
else
    [x, w, opp_ind] = deploy_all_octants(dim,x,w,adims);
end
w = (w /sum(w))*a_norm;
[nMtot, Sn, Kn] = compute_harmonics(dim, x, data.PnOrder, adims);
x = x(:,1:dim);
d2m = compute_d2m(nMtot, w, Sn);
m2d = compute_m2d(length(w), Sn, Kn, a_norm);
% Attach all matrices/vectors to data structures
% ----------------------------------------------
data.AngQuadNorm = a_norm;
data.NumberAngularDirections = length(w);
data.AngularDirections = x;
data.AngularWeights = w;
data.Opposite_Angular_Indices = opp_ind;
data.discrete_to_moment = d2m;
data.moment_to_discrete = m2d;
data.TotalFluxMoments = nMtot;
data.SphericalHarmonics = Sn;
data.MomentOrders = Kn;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_LDFE_quad(dim, level, adims)
% Quick Error Checks
% ------------------
if level < 0 || level > 7, error('LDFE level must be between 0 and 7.'); end
% Build Quadrature
% ----------------
octAngles = 4^(level+1);
numAngles = octAngles*2^dim;
angs = zeros(numAngles,3);
wts = zeros(numAngles,1);
lineno = 0;
for i=0:level-1
    lineno = lineno + 4^(i+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs, wts] = get_GLC_quad(dim, level, adims)
% Quick Error Checks
% ------------------

% Build Quadrature
% ----------------
ndir = (level*(level+2)/8)*2^dim;
az_nodes = zeros(level*(level+2)/8,1);
az_wt = zeros(level*(level+2)/8,1);
pol_wt = zeros(level/2,1);
cos_theta = zeros(level/2,1);
sin_theta = zeros(level/2,1);
[x,w] = lgwt(level,-1,1);
% Allocate output memory
angs = zeros(ndir,3);
wts = zeros(ndir,1);
for i=level/2:level-1
    cos_theta(i-level/2+1) = x(i);
    pol_wt(i-level/2+1) = w(i);
%     sin_theta(i-level/2+1) = sqrt(1-cos_theta(i-level/2+1)^2);
end
sin_theta = sqrt(1-cos_theta.^2);
pos = 1;
for i=0:level/2-1
    jmax = level/2-i;
    for j=0:level/2-i-1
        az_nodes(pos) = (pi/2)*j/jmax + (pi/4)/jmax;
        az_wt(pos) = (pi/2)/jmax;
        pos = pos + 1;
    end
end
pos = 1;offset = 0;
for i=0:level/2-1
    for j=0:level/2-i-1
        angs(pos,adims(1)) = sin_theta(i+1)*cos(az_nodes(j+1+offset));
        angs(pos,adims(2)) = sin_theta(i+1)*sin(az_nodes(j+1+offset));
        angs(pos,adims(3)) = cos_theta(i+1);
%         angs(pos,1) = sin_theta(i+1)*cos(az_nodes(j+1+offset));
%         angs(pos,2) = sin_theta(i+1)*sin(az_nodes(j+1+offset));
%         angs(pos,3) = cos_theta(i+1);
        wts(pos) = pol_wt(i+1)*az_wt(j+1+offset);
        pos = pos + 1;
    end
    offset = offset + level/2-i;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_PGLC_quad(dim, npolar, nazimuth, adims)
% Quick Error Checks
% ------------------
if npolar < 1, error('PGLC polar count must be >= 1.'); end
if nazimuth < 1, error('PGLC azimuthal count must be >= 1.'); end
if dim == 1, error('PGLC requires dimension > 1.'); end
% Build Quadrature
% ----------------
numAngles = npolar*nazimuth*2^dim;
nRoots = 2*npolar;
[x,w] = lgwt(nRoots,-1,1);
% Allocate output memory
angs = zeros(numAngles,3);
wts = zeros(numAngles,1);
% Get other local variables
delta_phi = pi/(4*nazimuth);
azim_weight = 2*delta_phi;
azim_angle = delta_phi;
for i=0:nazimuth-1
    for j=npolar:nRoots-1
        jj = j + 1;
        costheta = x(jj);
        sintheta = sqrt(1-costheta^2);
        iord = jj + (i-1)*npolar;
        angs(iord,adims) = [cos(azim_angle)*sintheta,sin(azim_angle)*sintheta,costheta];
        wts(iord) = azim_weight*w(jj);
    end
    azim_angle = azim_angle + 2*delta_phi;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_LS_quad(dim, level, adims)
% Quick Error Checks
% ------------------
if level > 24
    error('Error: Cannot go high thatn S24 for level-symmetric. Yields negative weights.')
end
% Build Quadrature
% ----------------
[angs,wts] = get_LS_local_quad(dim,level);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts] = get_TriGLC_quad(dim, npolar, nazimuth, adims)
% Quick Error Checks
% ------------------
if npolar < 1, error('TriPGLC polar count must be >= 1.'); end
if nazimuth < 1, error('TriPGLC azimuthal count must be >= 1.'); end
if dim == 1, error('TriPGLC requires dimension > 1.'); end
% Build Quadrature
% ----------------
octAngles = 0;
nnazimuth = nazimuth;
for i=1:npolar
    octAngles = octAngles + nnazimuth;
    nnazimuth = nnazimuth - 1;
end
numAngles = octAngles*2^dim;
nRoots = 2*npolar;
[x,w] = lgwt(nRoots,-1,1);
% Allocate memory
angs = zeros(numAngles,3);
wts  = zeros(numAngles,1);
% Loop through polar directions
dir = 1;
for j=0:npolar-1
    jj = j + 1;
    nazimuth_ = nazimuth;
    delta_phi = pi/(4.0*nazimuth_);
    azim_angle = delta_phi;
    azim_weight = 2*delta_phi;
    
    k = jj + npolar;
    costheta = x(k);
    sintheta = sqrt(1-costheta^2);
    for i=0:nazimuth-1
        angs(dir,adims) = [cos(azim_angle)*sintheta,sin(azim_angle)*sintheta,costheta];
%         angs(dir,:) = [cos(azim_angle)*sintheta,sin(azim_angle)*sintheta,costheta];
        wts(dir) = azim_weight*w(k);
        azim_angle = azim_angle + 2*delta_phi;
        dir = dir + 1;
    end
    nazimuth = nazimuth - 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [angs,wts,opp_ind] = deploy_all_octants(dim, angs, wts, adims)
numAngles = length(wts);
octAngles = numAngles / 2^dim;
opp_ind = zeros(numAngles,1);
if dim == 1
    for n=1:octAngles
        opp_ind(n) = n + 2*octAngles;
        opp_ind(n + 2*octAngles) = n;
    end
end
if dim == 2 || dim == 3
    for octant=2:4
        for n=1:octAngles
            m = (octant - 1)*octAngles + n;
            switch (octant)
                case(2)
                    angs(m,1) = -angs(n,1);
                    angs(m,2) =  angs(n,2);
                    angs(m,3) =  angs(n,3);
                case(3)
                    angs(m,1) = -angs(n,1);
                    angs(m,2) = -angs(n,2);
                    angs(m,3) =  angs(n,3);
                case(4)
                    angs(m,1) =  angs(n,1);
                    angs(m,2) = -angs(n,2);
                    angs(m,3) =  angs(n,3);
            end
            wts(m) = wts(n);
        end
    end
end
if dim == 3
    for n=1:numAngles/2
        angs(n+numAngles/2,1) =  angs(n,1);
        angs(n+numAngles/2,2) =  angs(n,2);
        angs(n+numAngles/2,3) = -angs(n,3);
        wts(n+numAngles/2)    =  wts(n);
    end
end
if dim == 2
    for n=0:octAngles-1
        nn = n + 1;
        opp_ind(nn) = nn+2*octAngles;
        opp_ind(nn+2*octAngles) = nn;
        opp_ind(nn+octAngles)   = nn+3*octAngles;
        opp_ind(nn+3*octAngles) = nn+octAngles;
    end
elseif dim == 3
    for n=0:octAngles-1
        nn = n + 1;
        opp_ind(nn) = nn+6*octAngles;
        opp_ind(nn+6*octAngles) = nn;
        opp_ind(nn+octAngles)   = nn+7*octAngles;
        opp_ind(nn+7*octAngles) = nn+octAngles;
        opp_ind(nn+2*octAngles) = nn+4*octAngles;
        opp_ind(nn+4*octAngles) = nn+2*octAngles;
        opp_ind(nn+3*octAngles) = nn+5*octAngles;
        opp_ind(nn+5*octAngles) = nn+3*octAngles;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opp_ind = determine_opposite_angle(angs)
na = size(angs,1);
opp_ind = zeros(na,1);
for m=1:na
    tang = angs(m,:);
    for mm=1:na
        if norm(tang - angs(mm,:)) < 1e-13
            opp_ind(m) = mm;
            break
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nMtot, Sn, Kn] = compute_harmonics(dim, angs, fMom, adims)
if dim == 1
    [nMtot, Sn, Kn] = compute_1D_harmonics(angs, fMom);
    return
end
% Determine Total Number of Flux Moments
% and determine Legendre coefficients.
c = 0;
for k=0:fMom
    switch(dim)
        case(1)
            lN = 0; uN = 0; step = 1;
        case(2)
            lN = -k; uN = k; step = 2;
        case(3)
            lN = -k; uN = k; step = 1;
    end
    for n=lN:step:uN
        c = c + 1;
        Kn(c,1) = k;
        Kn(c,2) = n;
    end
end
nMtot = size(Kn,1);
nA = size(angs,1);
Sn = zeros(nMtot, nA);
for a=1:nA
    L = evaluate_legendre(angs(a,adims(3)),fMom);
    for m=1:nMtot
        k = Kn(m,1); kk = k + 1;
        n = Kn(m,2); nn = n + 1;
        phi = atan(angs(a,adims(2)) / angs(a,adims(1)));
        if angs(a,adims(1)) < 0
            phi = phi + pi;
        end
        if n >= 0
            tf = cos(n*phi);
            pkn = (-1)^n * L(kk,nn);
        else
            tf = sin(n*phi);
            pkn = (-1)^(-n)*factorial(k+n)/factorial(k-n)*L(kk,-nn+1);
        end
        Sn(m,a) = sqrt((2-delta(n,0))*factorial(k-n)/factorial(k+n))*tf*pkn;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nMtot, Sn, Kn] = compute_1D_harmonics(angs, fMom)
nA = length(angs);
if fMom == 0
    nMtot = 1;
    Kn = [0, 0];
    Sn = ones(1, nA);
else
    nMtot = fMom+1;
    Kn = [(0:fMom)',(0:fMom)'];
    Sn = zeros(nMtot, nA);
    for m=0:fMom
        for q=1:nA
            Sn(m+1,q) = legendreP(m,angs(q));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d2m = compute_d2m(nM, w, Sn)
nA = length(w);
d2m = zeros(nM,nA);
for m=1:nM
    for a=1:nA
        d2m(m,a) = w(a)*Sn(m,a);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m2d = compute_m2d(nA, Sn, Kn, an)
nM = size(Kn,1);
m2d = zeros(nM,nA);
for m=1:nM
    for a=1:nA
        m2d(m,a) = (2*Kn(m,1) + 1) * Sn(m,a) / an;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = evaluate_legendre(x, order)
L = zeros(order+1);

for n=0:order
    nn = n + 1 ;
    for k=n:order
        kk = k + 1;
        if k==0
            L(kk,nn) = 1;
        else
            if k==n
                L(kk,nn) = (-1)^k*(1-x^2)^(k/2)*double_factorial(2*k-1);
            elseif k==n+1
                L(kk,nn) = x*(2*n+1)*L(kk-1,nn);
            elseif k>n+1
                L(kk,nn) = (x*(2*k-1)*L(kk-1,nn) - (k+n-1)*L(kk-2,nn))/(k-n);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = double_factorial(i)
if i == 0
    out = 1;
else
    out = 1;
    for j=i:-2:1
        out = out * j;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = delta(r, offset)
if r == offset
    out = 1;
else
    out = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
