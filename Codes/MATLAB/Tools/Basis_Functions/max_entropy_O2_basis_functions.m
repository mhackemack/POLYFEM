%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Maximum Entropy (2nd order) Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the basis function and gradient
%                   values using the Maximum Entropy Methodology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = max_entropy_O2_basis_functions( varargin )
nout = nargout;
grad_bool = false;
% Collect Input Arguments
% -----------------------
verts = varargin{1};
qx = varargin{2};
faces = varargin{3};
order = varargin{4};
nverts = varargin{5};
% Prepare Vertices and Dimensional Space
% --------------------------------------
dim = size(verts, 2);
nqx = size(qx, 1);
% Quick Error Checking
% --------------------
if order ~= 2, error('Should only be 2nd order in this functor.'); end
if dim ~= 2, error('Should only be 2D.'); end
% Allocate Matrix Memory
% ----------------------
if nout > 1, grad_bool = true; end
ntot = get_num_serendipity_points( dim, nverts, length(faces), order);
bout = zeros(nqx, ntot); pm_full = zeros(nqx, ntot); pp_full = zeros(nqx, ntot);
% Get Problem Preliminaries
% -------------------------
tol = 1e-13;
h = get_max_diamter( verts ); h0 = eye(dim)/h;
scaled_verts = (h0*verts')'; qx = (h0*qx')';
rva = get_vertex_differences( scaled_verts, qx );
[ca, dca, qa] = get_constraint_matrix(scaled_verts, rva, qx);
[ef, efder] = get_face_functions( scaled_verts, nverts, qx, dim, faces );

% rva = get_vertex_differences( verts, qx );
% [ca, dca, qa] = get_constraint_matrix(verts, rva, qx);
% [ef, efder] = get_face_functions( verts, nverts, qx, dim, faces );
[pf, pfder] = get_prior_functions( rva, dim, ef, efder, nverts );
dpx = get_number_lagrange_points( dim, order ); dpxz = zeros(dpx, 1);
% Perform all Newton Iterations
% -----------------------------
iters = zeros(nqx, 1);
lam = zeros(nqx, dpx);

for q=1:nqx
    xx = dpxz; converged = false;
    % Perform Newton Iterations
    while ~converged && iters(q) < 200
        % Check convergence
        [pm, pp] = get_entropy_functions( xx, pf, ca, q );
        g = get_F( pm, pp, ca, qa, q );
        norm_g = norm(g);
        if norm_g < tol; break; end
        % Get Hessian
        H = get_J( pm, pp, ca, q );
        % Get search direction
        dx = -H\g;
        % Get step length
        a = get_step_size( norm_g, xx, dx, pf, ca, q );
        % Update lambda multipliers
        xx = xx + a*dx;
        % Update Counters
        iters(q) = iters(q) + 1;
    end
    % Finalize Basis Function Values for this Quadrature Point
    % --------------------------------------------------------
    lam(q,:) = xx';
    [pm, pp] = get_entropy_functions( xx, pf, ca, q );
    pm_full(q,:) = pm; pp_full(q,:) = pp;
    bout(q,:) = pp - pm; %bout(q,bout(q,:) < 1e-13) = 0;
%     pm_full(q,:) = pm_full(q,:) / sum(bout(q,:));
%     pp_full(q,:) = pp_full(q,:) / sum(bout(q,:));
    bout(q,:) = bout(q,:) / sum(bout(q,:));
end

% Create Basis Function Gradients
% -------------------------------
if grad_bool
    gout = get_basis_gradients( lam, pf, pfder, pm_full, pp_full, ca, dca );
    for q=1:nqx
        gout(:,:,q) = gout(:,:,q)*h0;
    end
end
% Set Output Arguments
% --------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_max = get_max_diamter( verts )
nv = size(verts,1);
out_max = 0;
for i=1:nv
    vi = verts(i,:);
    for j=1:nv
        if i==j, continue; end
        h = norm(verts(j,:) - vi);
        if h > out_max, out_max = h; end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv = get_vertex_differences( verts, qx )
[nv, dim] = size(verts);
nvones = ones(nv,1);
nqx = size(qx, 1);
rv = zeros(2*nv, dim, nqx);
for q=1:nqx
    rv(1:nv,:,q) = verts - nvones*qx(q,:);
    for i=1:nv
        an = nv + i;
        ii = [i,mod(i,nv)+1];
%         fc = mean(verts(ii,:));
        fc = (verts(ii(1),:) + verts(ii(2),:))/2;
        rv(an,:,q) = fc - qx(q,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ef, efder] = get_face_functions(verts, nverts, qx, dim, faces)
nqx = size(qx, 1);
rv = verts(1:nverts,:);
ef = zeros(nqx, nverts);
efder = zeros(nverts, dim, nqx);
for q=1:nqx
    drv = [rv(:,1) - qx(q,1), rv(:,2) - qx(q,2)];
    for i=1:nverts
        ii = [i, mod(i,nverts)+1];
        norm1 = norm(drv(ii(1),:));
        norm2 = norm(drv(ii(2),:));
        ef(q,i) =  norm1 + norm2 - norm(verts(ii(2),:) - verts(ii(1),:));
        efder(i,:,q) = -drv(ii(1),:)/norm1 - drv(ii(2),:)/norm2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ca, dca, qa] = get_constraint_matrix(verts, rva, qx)
p = 6;
nqx = size(qx, 1); nv = size(verts, 1);
qa = [1;zeros(p-1,1)];
ca = zeros(p,2*nv,nqx);
dca = zeros(p,2,2*nv,nqx);
% constant
ca(1,:,:) = 1;
% linear
ca(2,:,:) = squeeze(rva(:,1,:));
ca(3,:,:) = squeeze(rva(:,2,:));
% quadratic
for q=1:nqx
    trv = rva(:,:,q);
    for i=1:nv
        ii = mod(i,nv) + 1;
        xba = verts(ii,:) - verts(i,:);
        an = i + nv;
        ca(4,i,q) = trv(i,1)*trv(i,1);
        ca(5,i,q) = trv(i,2)*trv(i,2);
        ca(6,i,q) = trv(i,1)*trv(i,2);
        ca(4,an,q) = trv(an,1)*trv(an,1) - xba(1)*xba(1)/4;
        ca(5,an,q) = trv(an,2)*trv(an,2) - xba(2)*xba(2)/4;
        ca(6,an,q) = trv(an,1)*trv(an,2) - xba(1)*xba(2)/4;
    end
end
% gradients
dca(2,1,:,:) = -1; dca(3,2,:,:) = -1;
for q=1:nqx
    trv = rva(:,:,q);
    for i=1:nv
        an = i + nv;
        dca(4,1,i,q) = -2*trv(i,1); dca(4,1,an,q) = -2*trv(an,1);
        dca(5,2,i,q) = -2*trv(i,2); dca(5,2,an,q) = -2*trv(an,2);
        dca(6,1,i,q) = -trv(i,2);   dca(6,1,an,q) = -trv(an,2);
        dca(6,2,i,q) = -trv(i,1);   dca(6,2,an,q) = -trv(an,1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w, wder] = get_prior_functions(rva, dim, ef, efder, nverts)
nx = size(rva, 3); nt = nverts*2; ot = ones(1,nt);
w = zeros(nx, nt); wder = zeros(nt, dim, nx);
for q=1:nx
    for i=1:nverts
        tefder = efder(:,:,q);
        if i==1
            ii = [nverts,i];
        else
            ii = [i-1,i];
        end
        an = i + nverts;
        w(q,i) = 1/(ef(q,ii(1)) * ef(q,ii(2)));
        w(q,an) = 1/ef(q,ii(2));
    end
    sumw = sum(w(q,:));
    % form temp gradients
    temp_dw = zeros(nt, dim);
    for i=1:nverts
        ain = i + nverts;
        if i==1
            ii = [nverts,1];
        else
            ii = [i-1,i];
        end
        temp_dw(i,:)  = -w(q,i)^2 * ( ef(q,ii(1))*tefder(ii(2),:) + ef(q,ii(2))*tefder(ii(1),:) );
        temp_dw(ain,:) = -w(q,ain)^2 * tefder(i,:);
    end
    % Form true prior weight function gradients
%     twder = wder(:,:,q);
    twder = sumw * temp_dw;
    for i=1:nt
        twder(i,:) = twder(i,:) - w(q,i)*ot*temp_dw;
%         twder(i,:) = sumw * temp_dw(i,:);
%         for j=1:nt
%             twder(i,:) = twder(i,:) - w(q,i) * temp_dw(j,:);
%         end
    end
    w(q,:) = w(q,:) / sumw;
    wder(:,:,q) = twder / sumw^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_number_lagrange_points( dim, order )
if dim == 2
    out = (order+1)*(order+2)/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = get_Z( x, pf, ca, q )
[pm, pp] = get_entropy_functions( x, pf, ca, q );
Z = x(1) + sum(pm + pp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pm, pp] = get_entropy_functions( x, pf, ca, q )
% [nx, nv] = size(pf);
% pm = zeros(nx, nv); pp = zeros(nx, nv);
% for q=1:nx
%     cat = ca(:,:,q);
%     w = pf(q,:);
%     for j=1:nv
%         fm = -1 + x(q,:)*cat(:,j);
%         fp = -1 - x(q,:)*cat(:,j);
%         pm(q,j) = w(j)*exp(fm);
%         pp(q,j) = w(j)*exp(fp);
%     end
% end
% pm = zeros(1, nv); pp = zeros(1, nv);
% cat = ca(:,:,q); 
% for j=1:nv
%     fm = -1 + xt*cat(:,j);
%     fp = -1 - xt*cat(:,j);
%     pm(j) = w(j)*exp(fm);
%     pp(j) = w(j)*exp(fp);
% end
% fm = -1 + xt*cat;
% fp = -1 - xt*cat;
% pm = w.*exp(fm);
% pp = w.*exp(fp);
w = pf(q,:);
xt = x';
xtc = xt*ca(:,:,q);
pm = w.*exp(-1 + xtc);
pp = w.*exp(-1 - xtc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = get_F( pm, pp, ca, qa, q )
% Z = pp - pm;
F = qa - ca(:,:,q)*(pp - pm)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = get_J( pm, pp, ca, q )
% [nx, nt] = size(pf);
% nt = length(pm);
% J = zeros(6);
cat = ca(:,:,q);
% for i=1:nt
%     J = J + (pp(i) + pm(i))*(cat(:,i)*cat(:,i)');
% end
ot = ones(6,1);
J = (cat.*(ot*(pm+pp))) * cat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_step_size( normg, x0, dx, pf, ca, q )
a = 1;
if normg > 10^(-4)
    a3 = a;
    g1 = get_Z(x0,pf,ca,q);
    g3 = get_Z(x0+a3*dx,pf,ca,q);
    % Loop until minimum energy aquired
    while g3 > g1 && a3 > 0.5*1e-13
        a3 = 0.5*a3;
        x = x0 + a3*dx;
        g3 = get_Z(x,pf,ca,q);
    end
    % Final step size calc
    if a3 > 0.5*1e-13
        a2 = 0.5*a3;
        x = x0 + a2*dx;
        g2 = get_Z(x,pf,ca,q);
        h1 = (g2 - g1)/a2;
        h2 = (g3 - g2)/(a3-a2);
        h3 = (h2 - h1)/a3;
        if abs(h3) > 1e-13
            a1 = 0.5*(a2 - h1/h3);
        else
            a1 = a3;
        end
        x = x0 + a1*dx;
        g0 = get_Z(x,pf,ca,q);
        if g0 < g3
            a = a1;
        else
            a = a3;
        end
    else
        a = a3;
    end
end
out = a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xx = vec_to_mat( x, nx, dim )
xx = zeros(nx, dim);
for i=1:nx
    xx(i,:) = x((i-1)*dim+1:i*dim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = mat_to_vec( x, nx, ~ )
x = reshape(x,nx,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grads = get_basis_gradients( x, pf, pfder, pm, pp, ca, dca )
[nqx, nt] = size(pf); zt = [1,1];
grads = zeros(nt, 2, nqx);
% Assemble Gradients
Z = ((pp - pm)./pf)'; pmp = (pm + pp)';
[H, A] = get_derivative_mats( x, pf, pfder, pm, pp, ca, dca );
for q=1:nqx
    tdca = dca(:,:,:,q);
    tHA = H(:,:,q)\A(:,:,q);
    grads(:,:,q) = (Z(:,q)*zt).*pfder(:,:,q);
    grads(:,:,q) = grads(:,:,q) + pmp(:,q)*zt.*(ca(:,:,q)'*tHA);
    for i=1:nt
%         grads(i,:,q) = grads(i,:,q) + pmp(i,q) * (ca(:,i,q)'*tHA);
        grads(i,:,q) = grads(i,:,q) - pmp(i,q) * (x(q,:)*tdca(:,:,i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H, A] = get_derivative_mats( x, pf, pfder, pm, pp, ca, dca )
[nqx, nt] = size(pf);
% allocate memory
H = zeros(6,6,nqx); A = zeros(6,2,nqx);
Z = ((pp - pm)./pf); pmp = (pm + pp);
% loop through quadrature points
for q=1:nqx
    H(:,:,q) = get_J( pm(q,:), pp(q,:), ca, q );
    for i=1:nt
%         A(:,:,q) = A(:,:,q) - Z(q,i)*ca(:,i,q)*pfder(i,:,q);
%         A(:,:,q) = A(:,:,q) + pmp(q,i)*ca(:,i,q)*(x(q,:)*dca(:,:,i,q));
%         A(:,:,q) = A(:,:,q) - (pp(q,i) - pm(q,i))*dca(:,:,i,q);
        
        A(:,:,q) = A(:,:,q) - Z(q,i)*ca(:,i,q)*pfder(i,:,q)...
                            + pmp(q,i)*ca(:,i,q)*(x(q,:)*dca(:,:,i,q))...
                            - (pp(q,i) - pm(q,i))*dca(:,:,i,q);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function grads = get_basis_gradients( x, pf, pfder, pm, pp, ca, dca )
% [nqx, nt] = size(pf);
% grads = zeros(nt, 2, nqx);
% [H, A] = get_derivative_mats( x, pf, pfder, pm, pp, ca, dca );
% % assemble gradients here
% Z = (pp - pm)./pf;
% for q=1:nqx
%     xx = x(q,:);
%     cat = ca(:,:,q)';
%     Hinv = inv(H(:,:,q));
%     for i=1:nt
%         grads(i,:,q) = Z(q,i)*pfder(i,:,q) + (pp(q,i) + pm(q,i))*(cat(i,:)*Hinv*A(:,:,q) - xx*dca(:,:,i,q));
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [H, A] = get_derivative_mats( x, pf, pfder, pm, pp, ca, dca )
% [nx, nt] = size(pf);
% H = zeros(6,6,nx); A = zeros(6,2,nx);
% % Loop through quadrature points
% Z = (pp - pm)./pf; b = pp - pm;
% for q=1:nx
%     xx = x(q,:);
%     H(:,:,q) = get_J( pf, pm, pp, ca, q );
%     tca = ca(:,:,q); tdca = dca(:,:,:,q);
%     AA = A(:,:,q);
%     for i=1:nt
%         AA = AA - Z(q,i)*tca(:,i)*pfder(i,:,q);
%         AA = AA + (pp(q,i) + pm(q,i))*tca(:,i)*(xx*tdca(:,:,i));
%         AA = AA - b(q,i)*tdca(:,:,i);
%     end
%     A(:,:,q) = AA;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%