function varargout = barycentric_serendipity(v, x, faces, basis)
% Get Input/Output Information
% ----------------------------
[nv, dim] = size(v);
nx = size(x,1);
nout = nargout;
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
% Get Basis Set Values and Transformations
% ------------------------------------------------------------------------------
[quad_pairs, qp1, qp2, nq] = get_quad_pairings(nv);
diag_pairs = get_diagonal_pairings(nv, quad_pairs);
% Get basis values
bvals = basis(v, x, faces, 1, nv);
grad_bool = false;
q_vals = zeros(nx, nq);
for i=1:nx
    q_vals(i,:) = bvals(i,qp1).*bvals(i,qp2);
end
% Get basis gradients if necessary
if nout == 2
    [bvals, bgrads] = basis(v, x, faces, 1, nv);
    grad_bool = true;
    q_grads = zeros(nq, dim, nx);
    q_vals = bvals(:,qp1).*bvals(:,qp2);
    for i=1:nx
        q_grads(:,:,i) = bgrads(qp1,:,i).*bgrads(qp2,:,i);
    end
end
A = get_quad_pairing_transformation(nv, v, vind, diag_pairs);
% B = get_lagrange_transformation(nv);
% Perform Transformations
% ------------------------------------------------------------------------------
ser_vals = zeros(nx, 2*nv);
ser_grads = zeros(2*nv, dim, nx);
for i=1:nx
    ser_vals(i,:) = A*q_vals(i,:)';
end
for i=1:nx
%     ser_vals(i,:) = (B*(A*q_vals(i,:)'))';
%     if grad_bool
%         for d=1:dim
%             ser_grads(:,d,i) = (B*(A*q_grads(:,d,i)))';
%         end
%     end
end
% Set Outputs
% -----------
ser_vals(:,nv+1:end) = 2*ser_vals(:,nv+1:end);
varargout{1} = ser_vals;
% if grad_bool
%     % Loop through all serendipity points
%     for q=1:nx
%         for i=1:2*nv
%             tt = [0,0];
%             t1 = quad_pairs(i,1); t2 = quad_pairs(i,2);
%             tt = tt + blin(q,t1)*glin(t2,:,q);
%             tt = tt + blin(q,t2)*glin(t1,:,q);
%         end
%     end
%     varargout{2} = ser_grads;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_quad_pairing_transformation(nv, v, quad_pairs, diag_pairs)
nd = size(diag_pairs,1);
I = eye(2*nv);
out = [I,get_a_prime_matrix(nd, nv, v, quad_pairs, diag_pairs)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_lagrange_transformation(nv)
I = eye(nv);
a = diag(-ones(nv,1),0) + diag(-ones(nv-1,1),-1); a(1,nv) = -1;
out = [I,a;zeros(nv),4*I];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qps, qp1, qp2, nq] = get_quad_pairings(nv)
nq = nv*(nv+1)/2;
out = zeros(nq,2);
for i=1:nv
    ii = [i, mod(i,nv)+1];
    out(i,:) = i;
    out(nv+i,:) = ii;
end
% Set Diagonal Terms
% ------------------
if nv > 3
    d = 2*nv;
    for i=1:nv-2
        if i==1
            jj = nv-1;
        else
            jj = nv;
        end
        for j=i+2:jj
            d = d + 1;
            out(d,:) = [i,j];
        end
    end
%     ntop = floor(nv/2) + rem(nv,2);
%     for i=1:ntop
%         if i ~= ntop
%             for j=i-2:-1:1
%                 d = d + 1;
%                 out(d,:) = [i,j];
%             end
%         end
%         for j=i+2:nv
%             if i==1 && j==nv, continue; end
%             d = d + 1;
%             out(d,:) = [i,j];
%         end
%     end
end
qps = out;
qp1 = out(:,1)';
qp2 = out(:,2)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_diagonal_pairings(nv, qps)
if nv < 4, out = []; return; end
nq = nv*((nv+1)/2-2);
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
function out = get_a_prime_matrix(nd, nv, v, vind, diag_pairs)
out = zeros(2*nv, nd);
for i=1:nd
    cvec = zeros(6, 1);
    % get rotation information for each diagonal
    dp = diag_pairs(i,:);
    vind1 = vind(dp(1),:);
    vind2 = vind(dp(2),:);
    vm = mean(v(dp(1:2),:));
    len = norm(diff(v(dp(1:2),:)))/2;
    vv = [v(:,1)-vm(1),v(:,2)-vm(2)]; vvv = vv'; 
    vvt = vv(dp(1),:);
    thet = acos(vv(dp(2),1)/len);
    tmat = [cos(thet), -sin(thet);  sin(thet), cos(thet)];
    vvt = tmat*vvt';
    if abs(vvt(2)) > 1e-14
        thet = 2*pi - thet;
        tmat = [cos(thet), -sin(thet);  sin(thet), cos(thet)];
    end
    for j=1:nv
        vv(j,:) = (tmat*vvv(:,j))';
    end
    % determine submatrix coefficients
    da = (vv(vind1(1),1)*vv(vind1(2),2) - vv(vind1(1),2)*vv(vind1(2),1)) / (vv(vind1(1),2) - vv(vind1(2),2)) / len;
    db = (vv(vind2(2),1)*vv(vind2(1),2) - vv(vind2(2),2)*vv(vind2(1),1)) / (vv(vind2(1),2) - vv(vind2(2),2)) / len;
    s = 2 / (2 - (da + db));
    
    cvec(1) = -s*(1+da);
    cvec(2) = -s*(1+db);
    
    A34 = [1, 1; vv(vind1(1),2), vv(vind1(2),2)];
    A56 = [1, 1; vv(vind2(1),2), vv(vind2(2),2)];
    b34 = s*[1;da*vv(dp(1),2)];
    b56 = s*[1;db*vv(dp(2),2)];
    cvec(3:4) = A34\b34;
    cvec(5:6) = A56\b56;
    % set into global submatrix
    out(dp,i) = cvec;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%