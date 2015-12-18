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
%   A Matrix Notes: To determine the quadratic to serendipity transformation
%                   matrix (xi=A*mu), we utilize the Moore-Penrose psuedoinverse.
%                   Eliminating each diagonal entry provides an underdetermined
%                   system of equations. Given a (6x2v) system of constraint
%                   equations, L, the pseudoinverse is L'*(L*L')^(-1).
%
%                   https://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse
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
vind = get_adjacent_vertices(nverts);
[ser_verts, ser_nodes] = get_serendipity_nodes(nverts, scaled_verts, faces);
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
% A = get_quad_pairing_transformation(nverts, scaled_verts, vind, diag_pairs);
A = get_quad_pairing_transformation(ser_verts, ser_nodes, quad_pairs, diag_pairs, vind);
q_vals = blin(:,quad_pairs(:,1)).*blin(:,quad_pairs(:,2));
for q=1:nqx
    bout(q,:) = A*q_vals(q,:)';
end
% Create Basis Function Gradients
% ------------------------------------------------------------------------------
if grad_bool
    gfact = [ones(1,nverts),2*ones(1,nverts)];
    % Loop through points and calculate gradient values
    for q=1:nqx
        for i=1:ntot
            tt = [0,0];
            t1 = quad_pairs(i,1); t2 = quad_pairs(i,2);
            tt = tt + blin(q,t1)*glin(t2,:,q);
            tt = tt + blin(q,t2)*glin(t1,:,q);
            for j=1:num_dp
                t1 = quad_pairs(ntot+j,1); t2 = quad_pairs(ntot+j,2);
                aa = A(i,ntot+j);
                tt = tt + aa*blin(q,t1)*glin(t2,:,q);
                tt = tt + aa*blin(q,t2)*glin(t1,:,q);
            end
            gout(i,:,q) = gfact(i)*tt;
        end
        gout(:,:,q) = gout(:,:,q)*h0;
    end
end
% Set Output Arguments
% ------------------------------------------------------------------------------
bout(:,nverts+1:end) = 2*bout(:,nverts+1:end);
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
function vind = get_adjacent_vertices(nv)
vind = zeros(nv,2);
for i=1:nv
    if i == 1
        ii = [nv,i+1];
    elseif i == nv
        ii = [i-1,1];
    else
        ii = [i-1, i+1];
    end
    vind(i,:) = ii;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ser_verts, ser_nodes] = get_serendipity_nodes(nverts, verts, faces)
dim = size(verts,2);
ser_verts = zeros(2*nverts, dim);
ser_nodes = zeros(2*nverts, 2); ser_nodes(1:nverts,:) = [(1:nverts)',(1:nverts)'];
ser_verts(1:nverts,:) = verts;
for f=1:length(faces)
    fv = faces{f};
    ser_verts(nverts+f,:) = (verts(fv(1),:) + verts(fv(2),:))/2;
    ser_nodes(nverts+f,:) = fv;
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
% function A = get_quad_pairing_transformation(verts, ser_nodes, quad_pairs, diag_pairs, vind)
% nv = size(verts,1); nqp = size(diag_pairs,1);
% A = zeros(nv, size(quad_pairs,1)); A(1:nv,1:nv) = eye(nv);
% tL = zeros(6,nv); %tq = zeros(6,1);
% % tL = zeros(6); tq = zeros(6,1);
% % Loop through interior diagonal pairs
% d = nv;
% for i=1:nqp
%     d = d + 1;
%     va = verts(quad_pairs(d,1),:);
%     vb = verts(quad_pairs(d,2),:);
%     L = tL;
%     tdp = diag_pairs(i,:);
%     % Loop through all serendipity nodes
%     for j=1:6
%         v1 = verts(ser_nodes(j,1),:); v2 = verts(ser_nodes(j,2),:);
%         t = v1'*v2 + v2'*v1;
%         % c-constraint
%         L(1,j) = 2;
%         % x-constraint
%         L(2,j) = 2*verts(j,1);
%         % y-constraint
%         L(3,j) = 2*verts(j,2);
%         % x_x-constraint
%         L(4,j) = t(1,1);
%         % y_y-constraint
%         L(5,j) = t(2,2);
%         % x_y-constraint
%         L(6,j) = t(1,2);
%     end
%     % Apply right-hand side
%     q = [2;(va(1)+vb(1));(va(2)+vb(2));2*va(1)*vb(1);2*va(2)*vb(2);(va(1)*vb(2)+va(2)*vb(1))];
%     t = L'*((L*L')\q); t(abs(t) < 1e-14) = 0;
%     A(:,d) = t;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = get_quad_pairing_transformation(verts, ser_nodes, quad_pairs, diag_pairs, vind)
nv = size(verts,1); nqp = size(diag_pairs,1);
A = zeros(nv, nqp); A(1:nv,1:nv) = eye(nv);
tL = zeros(6); tq = zeros(6,1);
% Loop through interior diagonal pairs
d = nv;
for i=1:nqp
    d = d + 1;
    va = verts(quad_pairs(d,1),:);
    vb = verts(quad_pairs(d,2),:);
    tdp = diag_pairs(i,:);
    L = tL; q = tq;
    % c-constraint
    L(1,1:2) = 1; L(1,3:end) = 1; q(1) = 2;
    % x-constraint
    L(2,:) = verts(tdp,1); L(2,3:end) = 2*L(2,3:end);
    q(2) = (va(1) + vb(1));
    % y-constraint
    L(3,:) = verts(tdp,2); L(3,3:end) = 2*L(2,3:end);
    q(3) = (va(2) + vb(2));
    % x_x-constraint
    L(4,:) = (verts(ser_nodes(tdp,1),1).*verts(ser_nodes(tdp,2),1))';
    q(4) = va(1)*vb(1);
    % y_y-constraint
    L(5,:) = (verts(ser_nodes(tdp,1),2).*verts(ser_nodes(tdp,2),2))';
    q(5) = va(2)*vb(2);
    % x_y-constraint
    L(6,:) = ((verts(ser_nodes(tdp,1),1).*verts(ser_nodes(tdp,2),2))')/2 + ...
             ((verts(ser_nodes(tdp,1),2).*verts(ser_nodes(tdp,2),1))')/2;
    q(6) = (va(1)*vb(2) + va(2)*vb(1))/2;
    % Solve and apply constraint
    t = L'*((L*L')\q); t(abs(t) < 1e-14) = 0;
    A(tdp,d) = t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = get_quad_pairing_transformation(nv, v, quad_pairs, diag_pairs)
% nd = size(diag_pairs,1);
% I = eye(2*nv);
% out = [I,get_a_prime_matrix(nd, nv, v, quad_pairs, diag_pairs)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = get_a_prime_matrix(nd, nv, v, vind, diag_pairs)
% out = zeros(2*nv, nd);
% for i=1:nd
%     cvec = zeros(6, 1);
%     % get rotation information for each diagonal
%     dp = diag_pairs(i,:);
%     vind1 = vind(dp(1),:);
%     vind2 = vind(dp(2),:);
%     vm = mean(v(dp(1:2),:));
%     len = norm(diff(v(dp(1:2),:)))/2;
%     vv = [v(:,1)-vm(1),v(:,2)-vm(2)]; vvv = vv';
%     thet = acos(vv(dp(1),1)/len);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %     |  
%     %  4  |  3
%     %     |
%     %-----------
%     %     |
%     %  1  |  2
%     %     |
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if vv(dp(1),1) <= 0 && vv(dp(1),2) <= 0
% %         t = -(pi-thet);
%         thet = thet - pi;
%     elseif vv(dp(1),1) > 0 && vv(dp(1),2) <= 0
%         thet = thet + pi;
%     elseif vv(dp(1),1) > 0 && vv(dp(1),2) > 0
%         thet = pi - thet;
%     elseif vv(dp(1),1) <= 0 && vv(dp(1),2) > 0
%         thet = pi - thet;
%     end
%     tmat = [cos(thet), -sin(thet);  sin(thet), cos(thet)];
% %     if thet1 > pi/2
% %         tmat = [cos(thet1), -sin(thet1);  sin(thet1), cos(thet1)];
% %     else
% %         tmat = [cos(thet2), -sin(thet2);  sin(thet2), cos(thet2)];
% %     end
% %     vvt = tmat*vv(dp(1),:)';
% %     thet = acos(vv(dp(2),1)/len);
% %     tmat = [cos(thet), -sin(thet);  sin(thet), cos(thet)];
% %     vvt = tmat*vvt';
% %     if abs(vvt(2)) > 1e-14
% %         thet = 2*pi - thet;
% %         tmat = [cos(thet), -sin(thet);  sin(thet), cos(thet)];
% %     end
%     for j=1:nv
%         vv(j,:) = (tmat*vvv(:,j))';
%     end
%     % determine submatrix coefficients
%     da = (vv(vind1(1),1)*vv(vind1(2),2) - vv(vind1(1),2)*vv(vind1(2),1)) / (vv(vind1(1),2) - vv(vind1(2),2)) / len;
%     db = (vv(vind2(2),1)*vv(vind2(1),2) - vv(vind2(2),2)*vv(vind2(1),1)) / (vv(vind2(1),2) - vv(vind2(2),2)) / len;
%     s = 2 / (2 - (da + db));
%     
%     cvec(1) = -s*(1+da);
%     cvec(2) = -s*(1+db);
%     
%     A34 = [1, 1; vv(vind1(1),2), vv(vind1(2),2)];
%     A56 = [1, 1; vv(vind2(1),2), vv(vind2(2),2)];
%     b34 = s*[1;da*vv(dp(1),2)];
%     b56 = s*[1;db*vv(dp(2),2)];
%     cvec(3:4) = A34\b34;
%     cvec(5:6) = A56\b56;
%     % set into global submatrix
%     out(dp,i) = cvec;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%