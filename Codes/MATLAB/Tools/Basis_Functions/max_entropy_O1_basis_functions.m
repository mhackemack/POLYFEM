%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Maximum Entropy (1st order) Basis Function Generator
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
function varargout = max_entropy_O1_basis_functions( varargin )
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
if order ~= 1, error('Should only be 1st order in this functor.'); end
if dim ~= 2, error('We enforce only 2D MAXENT coordinates at this time.'); end
% Allocate Matrix Memory
% ----------------------
if nout > 1, grad_bool = true; end
ntot = get_num_serendipity_points( dim, nverts, length(faces), order);
bout = zeros(nqx, ntot);
% Get Problem Preliminaries
% -------------------------
tol = 1e-12;
h = get_max_diamter( verts ); h0 = eye(dim)/h;
scaled_verts = (h0*verts')'; qx = (h0*qx')';
h_n = get_min_nodal_spacing( scaled_verts );
rva = get_vertex_differences( scaled_verts, qx );
[ef, efder] = get_face_functions(scaled_verts, nverts, qx, dim, faces);
[pf, pfder] = get_prior_functions( h_n, rva, dim, ef, efder );
dpx = get_number_lagrange_points( dim, order ); dpxz = zeros(dpx, 1);
[inp, onp] = inpolygon(qx(:,1), qx(:,2), scaled_verts(:,1), scaled_verts(:,2));
% Perform all Newton Iterations
% -----------------------------
iters = zeros(nqx, 1);
% Loop through Quadrature Nodes - only evaluate interior nodes
for q=1:nqx
    if ~inp(q) || onp(q), continue; end
    xx = dpxz; converged = false;
    % Perform Newton Iterations
    while ~converged
        % Check convergence
        Z = get_Z( xx, rva, pf, q );
        g = get_F( Z, rva, pf, q );
        norm_g = norm(g);
        if norm_g < tol; break; end
        % Get Hessian
        H = get_J( Z, rva, pf, q );
        % Get search direction
        dx = -H\g;
        % Get step length
        a = get_step_size( norm_g, xx, dx, rva, pf, q );
        % Update lambda multipliers
        xx = xx + a*dx;
        % Update Counters
        iters(q) = iters(q) + 1;
    end
    % Finalize Basis Function Values for this Quadrature Point
    % --------------------------------------------------------
    Z = get_Z( xx, rva, pf,q )';
    bout(q,:) = Z/sum(Z);
end
% Evaluate boundary nodes
if sum(onp) > 0
    for f=1:length(faces)
        fv = faces{f};
        f1 = scaled_verts(fv(1),:); f2 = scaled_verts(fv(2),:);
        len = norm(f2 - f1);
        for q=1:nqx
            if onp(q)
                tqx = qx(q,:);
                d1 = norm(f1-tqx); d2 = norm(f2-tqx);
                if abs(d1+d2-len) < 1e-12
                    bout(q,fv) = [d2,d1]/len;
                    onp(q) = false;
                end
            end
        end
    end
end
% Create Basis Function Gradients
% -------------------------------
if grad_bool
    gout = get_basis_gradients( scaled_verts, pf, pfder, rva, bout, faces );
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
function out = get_min_nodal_spacing( verts )
[nv, dim] = size(verts);
out = zeros(nv, 1);
for i=1:nv
    d = sqrt(sum((ones(nv,1)*verts(i,:)-verts).^2,2)).'; d(i) = -1;
    d = sort(d);
    if dim == 1
        out(i) = verts(2) - verts(1);
    elseif dim == 2
        out(i) = d(3); 
    else
        out(i) = d(4);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv = get_vertex_differences( verts, qx )
[nv, dim] = size(verts);
nvones = ones(nv,1);
nqx = size(qx, 1);
rv = zeros(nv, dim, nqx);
for q=1:nqx
    rv(:,:,q) = verts - nvones*qx(q,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ef, efder] = get_face_functions(verts, nverts, qx, dim, faces)
nqx = size(qx, 1); nf = length(faces);
rv = verts(1:nverts,:);
ef = zeros(nqx, nf);
efder = zeros(nf, dim, nqx);
if dim == 2
    for q=1:nqx
        tefder = efder(:,:,q);
        drv = [rv(:,1) - qx(q,1), rv(:,2) - qx(q,2)];
        for i=1:nf
            ii = [i, mod(i,nverts)+1];
            norm1 = norm(drv(ii(1),:));
            norm2 = norm(drv(ii(2),:));
            ef(q,i) =  norm1 + norm2 - norm(diff(verts(ii,:)));
            tefder(i,:) = -drv(ii(1),:)/norm1 - drv(ii(2),:)/norm2;
        end
        efder(:,:,q) = tefder;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w, wder] = get_prior_functions(h_n, rva, dim, ef, efder)
nt = size(rva, 1); nqx = size(rva, 3);
w = zeros(nqx, nt); wder = zeros(nt, dim, nqx);
if dim == 2
    for q=1:nqx
        tefder = efder(:,:,q);
        for i=1:nt
            if i==1
                ii = [nt,1];
            else
                ii = [i-1,i];
            end
            w(q,i) = 1 / (ef(q,ii(1)) * ef(q,ii(2)));
        end
        sumw = sum(w(q,:));
        twder = wder(:,:,q);
        % Compute prior function gradients
        for i=1:nt
            if i==1
                ii = [nt,1];
            else
                ii = [i-1,i];
            end
            twder(i,:) = (-sumw/(ef(q,ii(1)) * ef(q,ii(2)))^2) * ( ef(q,ii(2))*tefder(ii(1),:) + ef(q,ii(1))*tefder(ii(2),:) );
            for j=1:nt
                if j==1
                    jj = [nt,1];
                else
                    jj = [j-1,j];
                end
                tder = -(1 / (ef(q,jj(1)) * ef(q,jj(2))))^2 * ( ef(q,jj(2))*tefder(jj(1),:) + ef(q,jj(1))*tefder(jj(2),:) );
                twder(i,:) = twder(i,:) - w(q,i)*tder;
            end
        end
        w(q,:) = w(q,:) / sumw;
        wder(:,:,q) = twder(:,:) / sumw^2;
    end
else
    g = 1.2; bt = g./(h_n.^2);
    for q=1:nqx
        trv = rva(:,:,q); 
        w(q,:) = exp(-bt.*(sum((trv.').^2).'))';
        wder(:,:,q)=2*((w(q,:)'*ones(1,dim)).*trv).*(bt*ones(1,dim)); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_number_lagrange_points( dim, order )
if order == 1
    out = dim;
elseif order == 2
    if dim == 2
        out = 6;
    elseif dim == 3
        out = 11;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = get_Z( x, rv, pf, q)
% [nx, nv] = size(pf);


% Z = zeros(nx, nv);
% for q=1:nx
%     trv = rv(:,:,q);
%     for j=1:nv
%         Z(q,j) = pf(q,j) * exp(-(x(q,:)*trv(j,:)'));
%     end
% end


% Z = zeros(1, nv);
% trv = rv(:,:,q)';
Z = pf(q,:).*exp(-(x'*rv(:,:,q)'));
% for j=1:nv
%     Z(j) = pf(q,j) * exp(-(trv(j,:)*x));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = get_logpartitionfunction( x, rv, pf, q)
Z = log(sum(get_Z( x, rv, pf, q)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = get_F( Z, rv, pf, q )
% ntot = length(dx);
% dim = size(rv,2);
% dpx = dim;
% [nx, nv] = size(pf);


% FF = zeros(nx, dim);
% for i=1:nx
%     FF(i,:) = -1.0*( Z(i,:) * rv(:,:,i)) / sum(Z(i,:));
% end
% % Set into full Gradient Vector
% F = zeros(ntot, 1);
% for i=1:nx
%     ii = (i-1)*dpx+1:i*dpx;
%     F(ii) = FF(i,:)';
% end



FF = -1.0*( Z * rv(:,:,q)) / sum(Z);
F =  FF';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = get_J( Z, rv, pf, q )
% ntot = length(dx);
dim = size(rv,2);
[nx, nv] = size(pf);


% vones = ones(nv,1);
% JJ = zeros(dim, dim, nx);
% for i=1:nx
%     trv = rv(:,:,i);
%     sumZ = Z(i,:)*vones;
%     JJJ = JJ(:,:,i);
%     for j=1:nv
%         ttrv = trv(j,:)'*trv(j,:);
%         JJJ = JJJ + (Z(i,j)/sumZ)*ttrv;
%     end
%     JJ(:,:,i) = JJJ;
% end
% % Set into full Hessian Matrix
% J = zeros(ntot, ntot);
% for i=1:nx
%     ii = (i-1)*dpx+1:i*dpx;
%     J(ii,ii) = JJ(:,:,i);
% end



J = zeros(dim);
trv = rv(:,:,q);
sumZ = sum(Z);
for j=1:nv
    ttrv = trv(j,:)'*trv(j,:);
    J = J + (Z(j)/sumZ)*ttrv;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_step_size( normg, x0, dx, rv, pf, q )
a = 1;
if normg > 10^(-4)
    a3 = a;
    g1 = get_logpartitionfunction(x0,rv,pf,q);
    g3 = get_logpartitionfunction(x0+a3*dx,rv,pf,q);
    % Loop until minimum energy aquired
    while g3 > g1 && a3 > 0.5*1e-13
        a3 = 0.5*a3;
        x = x0 + a3*dx;
        g3 = get_logpartitionfunction(x,rv,pf,q);
    end
    % Final step size calc
    if a3 > 0.5*1e-13
        a2 = 0.5*a3;
        x = x0 + a2*dx;
        g2 = get_logpartitionfunction(x,rv,pf,q);
        h1 = (g2 - g1)/a2;
        h2 = (g3 - g2)/(a3-a2);
        h3 = (h2 - h1)/a3;
        if abs(h3) > 1e-13
            a1 = 0.5*(a2 - h1/h3);
        else
            a1 = a3;
        end
        x = x0 + a1*dx;
        g0 = get_logpartitionfunction(x,rv,pf,q);
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
function grads = get_basis_gradients( verts, pf, pfder, rva, b, faces )
[nv, dim] = size(verts);
nx = size(b, 1);
grads = zeros(nv, dim, nx);
[H,MA,A,MC] = get_derivative_mats( pf, pfder, rva, b, nv, dim, nx);
for i=1:nx
    phi = b(i,:)';
    phid = phi*ones(1,dim);
    invH = inv(H(:,:,i));
    grads(:,:,i) = phid.*(rva(:,:,i)*(invH - invH*A(:,:,i))) + phid.*MA(:,:,i) - phi*MC(i,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,MA,A,MC] = get_derivative_mats( pf, pfder, rva, b, nv, dim, nx)
H = zeros(dim, dim, nx);
MA = zeros(nv, dim, nx);
A = zeros(dim, dim, nx);
MC = zeros(nx, dim);
for i=1:nx
    trv = rva(:,:,i);
    phi = b(i,:);
    MA(:,:,i) = (ones(nv,1)*(1./pf(i,:))) * pfder(:,:,i);
    for d=1:dim
        MC(i,d) = phi*MA(:,d,i);
        for dd=1:dim
            H(d,dd,i) = phi* (trv(:,d).*trv(:,dd));
            A(d,dd,i) = (phi'.*MA(:,dd))'*trv(:,d);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%