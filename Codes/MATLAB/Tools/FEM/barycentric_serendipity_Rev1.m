%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Quadratic Serendipity Generator (Rev1)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input Space:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = barycentric_serendipity_Rev1(varargin)
nout = nargout;
grad_bool = false;
% Collect Input Arguments
% ------------------------------------------------------------------------------
verts = varargin{1};
qx = varargin{2};
faces = varargin{3};
basis = varargin{4};
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[nverts,dim] = size(verts); ntot = 2*nverts;
nqx = size(qx, 1);
% Quick Error Checking
% ------------------------------------------------------------------------------
if dim ~= 2, error('Should only be 2D.'); end
% Allocate Matrix Memory
% ------------------------------------------------------------------------------
bout = zeros(nqx, ntot); gout = zeros(ntot,dim,nqx);
% Get Problem Preliminaries
% ------------------------------------------------------------------------------
if nout > 1, grad_bool = true; end
h = get_max_diamter( verts ); h0 = eye(dim)/h;
scaled_verts = (h0*verts')'; qx = (h0*qx')';
% rva = get_vertex_differences( scaled_verts, qx );
quad_pairs = get_quad_pairings(nverts);
diag_pairs = get_diag_pairings(nverts, quad_pairs); num_dp = size(diag_pairs,1);
% Retrieve Linear Basis Functions
% ------------------------------------------------------------------------------
if grad_bool
    [blin, glin] = basis(scaled_verts,qx,faces,1,nverts);
else
    blin = basis(scaled_verts,qx,faces,1,nverts);
end
% Build Quadratic Serendipity Basis Function Space
% ------------------------------------------------------------------------------
A = get_quad_pairing_transformation(scaled_verts, quad_pairs, diag_pairs);
q_vals = blin(:,quad_pairs(:,1)).*blin(:,quad_pairs(:,2));
for q=1:nqx
    bout(q,:) = A*q_vals(q,:)';
end
% Create Basis Function Gradients
% ------------------------------------------------------------------------------
if grad_bool
    % Loop through points and calculate gradient values
    for q=1:nqx
        for i=1:ntot
            tt = [0,0];
            t1 = quad_pairs(i,1); t2 = quad_pairs(i,2);
            tt = tt + blin(q,t1)*glin(t2,:,q);
            tt = tt + blin(q,t2)*glin(t1,:,q);
            for j=1:num_dp
                tdp = diag_pairs(j,:);
                for k=1:length(tdp)
                    t1 = quad_pairs(tdp(k),1); t2 = quad_pairs(tdp(k),2);
                    aa = A(tdp(k),ntot+j);
                    tt = tt + aa*blin(q,t1)*glin(t2,:,q);
                    tt = tt + aa*blin(q,t2)*glin(t1,:,q);
                end
            end
            gout(i,:,q) = tt;
            gout(:,:,q) = gout(:,:,q)*h0;
        end
    end
end
% Set Output Arguments
% ------------------------------------------------------------------------------
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
        fc = (verts(ii(1),:) + verts(ii(2),:))/2;
        rv(an,:,q) = fc - qx(q,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quad_pairs = get_quad_pairings(nv)
ntot = 2*nv+nv*(nv-3)/2;
quad_pairs = zeros(ntot,2);
% Loop through vertices and assign V/E terms
d = 2*nv;
for i=1:nv
    quad_pairs(i,:) = i;
    quad_pairs(i+nv,:) = [i,mod(i,nv)+1];
    if i>nv-2, continue; end
    if i==1
        jj = nv-1;
    else
        jj = nv;
    end
    for j=i+2:jj
        d = d + 1;
        quad_pairs(d,:) = [i,j];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_diag_pairings(nv, qps)
nq = nv*(nv-3)/2;
out = zeros(nq,6);
d = 2*nv;
for i=1:nq
    d = d + 1;
    tqp = qps(d,:);
    out(i,1:2) = tqp;
    % 1st vertex
    if tqp(1) == 1
        out(i,3) = 2*nv;
    else
        out(i,3) = nv + tqp(1) - 1;
    end
    out(i,4) = nv + tqp(1);
    % second vertex
    if tqp(2) == 1
        out(i,5) = 2*nv;
    else
        out(i,5) = nv + tqp(2) - 1;
    end
    out(i,6) = nv + tqp(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = get_quad_pairing_transformation(verts, quad_pairs, diag_pairs)
nv = size(verts,1); nqp = size(quad_pairs,1);
A = zeros(2*nv, nqp); A(1:2*nv,1:2*nv) = eye(2*nv);
tL = zeros(6); tq = zeros(6,1);
% Loop through interior diagonal pairs
d = 2*nv;
for i=1:nqp
    d = d + 1;
    va = verts(quad_pairs(d,1),:);
    vb = verts(quad_pairs(d,2),:);
    tdp = diag_pairs(i,:);
    L = tL; q = tq;
    % c-constraint
    L(1,:) = 1; q(1) = 1;
    % x-constraint
%     L(2,:) = ; q(2) = (va(1) + vb(1))/2;
    % y-constraint
%     L(3,:) = ; q(3) = (va(2) + vb(2))/2;
    % x_x-constraint
    q(4) = va(1)*vb(1);
    % y_y-constraint
    q(5) = va(2)*vb(2);
    % x_y-constraint
    q(6) = (va(1)*vb(2) + va(2)*vb(1))/2;
    % Solve and apply constraint
    A(tdp,d) = L\q;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%