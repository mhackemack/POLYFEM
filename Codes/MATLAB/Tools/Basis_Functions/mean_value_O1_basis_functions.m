%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Linear Mean Value Basis Function Generator
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
function varargout = mean_value_O1_basis_functions(varargin)
nout = nargout;
grad_bool = false; if nout > 1, grad_bool = true; end
% Collect Input Arguments
% ------------------------------------------------------------------------------
verts = varargin{1};
qx = varargin{2};
faces = varargin{3};
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[nverts, dim] = size(verts);
nqx = size(qx, 1);
% Allocate Matrix Memory
% ------------------------------------------------------------------------------
if nout > 1, grad_bool = true; end
ntot = get_num_serendipity_points( dim, nverts, length(faces), 1);
bout = zeros(nqx, ntot);
if grad_bool, gout = zeros(ntot, dim, nqx); end
% Get Problem Preliminaries
% ------------------------------------------------------------------------------
[rv, sv] = get_vertex_differences(verts, qx);
% Build Basis Function Sets
% ------------------------------------------------------------------------------
if dim == 2
    [bout, R] = mv2D(bout, rv, sv, grad_bool);
elseif dim == 3
    [bout, R] = mv3D(bout, rv, sv, grad_bool);
end
if grad_bool, gout = get_basis_gradients(bout, R); end
% Set Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rv, sv] = get_vertex_differences( verts, qx )
[nv, dim] = size(verts);
nqx = size(qx, 1);
rv = zeros(nv, dim, nqx);
sv = zeros(nv, nqx);
zn = ones(nv, 1);
zones = ones(dim,1);
for q=1:nqx
    rv(:,:,q) = verts - zn*qx(q,:);
    rrv = rv(:,:,q).*rv(:,:,q);
    sv(:,q) = sqrt(rrv*zones);
end
sv = sv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bout, R] = mv2D(bout, rv, sv, grad_bool)
[nqx, nverts] = size(sv);
R = zeros(nverts, 2, nqx);
zn = zeros(1, nverts); ont = ones(nverts, 1);
A = zeros(nqx, nverts);
D = zeros(nqx, nverts);
v_bool = false(nqx,1);
e_bool = false(nqx,1); e_num = zeros(nqx,1);
% Check on vertex or face
for q=1:nqx
    rrv = rv(:,:,q);
    rrvt = rrv';
    for i=1:nverts
        % Check if on vertex - gradient is undefined at this point
        if abs(sv(q,i)) < 1e-13
            bout(q,i) = 1.0;
            v_bool(q) = true;
            break;
        end
        % Check if on edge
        ii = mod(i,nverts) + 1;
        J = [rrv(i,:);rrv(ii,:)];
        A(q,i) = (J(1,1)*J(2,2)-J(2,1)*J(1,2))/2;
        D(q,i) = rrv(i,:)*rrvt(:,ii);
        if abs(A(q,i)) < 1e-13 && D(q,i) < 0
            rr = sv(q,i) + sv(q,ii);
            bout(q,i) = sv(q,ii)/rr;
            bout(q,ii) = sv(q,i)/rr;
            e_bool(q) = true;
            e_num(q) = i;
            break;
        end
    end
    % Compute edge gradient ratio
    if grad_bool && e_bool(q)
        ee = [e_num(q),mod(e_num(q),nverts)+1];
        % Loop through all vertices
        for i=1:nverts
            im = mod(i-2,nverts) + 1;
            ip = mod(i,nverts) + 1;
            em = rrv(im,:)/sv(q,im);
            ei = rrv(i,:)/sv(q,i);
            ep = rrv(ip,:)/sv(q,ip);
            % First edge vertex
            if i == ee(1)
                R(i,:,q) = ei/sv(q,i);
            % Second edge vertex
            elseif i == ee(2)
                R(i,:,q) = ei/sv(q,i);
            % Non-adjacent vertices
            else
                cm = em/sv(q,im) - ei/sv(q,i);
                ci = ei/sv(q,i) - ep/sv(q,ip);
                am = acos((rrv(im,:)*rrvt(:,i)) / (sv(q,im)*sv(q,i)));
                ai = acos((rrv(ip,:)*rrvt(:,i)) / (sv(q,ip)*sv(q,i)));
                tm = tan(am/2);
                ti = tan(ai/2);
                R(i,:,q) = ei/sv(q,i);
                if abs(sin(am)) > 1e-13
                    R(i,:,q) = R(i,:,q) + tm/(tm+ti)*[-cm(2),cm(1)]/sin(am);
                end
                if abs(sin(ai)) > 1e-13
                    R(i,:,q) = R(i,:,q) + ti/(tm+ti)*[-ci(2),ci(1)]/sin(ai);
                end
            end
        end
    end
end
% Compute volume coordinates
for q=1:nqx
    if v_bool(q) || e_bool(q), continue; end
    rrv = rv(:,:,q);
    rrvt = rrv';
    w = zn;
    for i=1:nverts
        im = mod(i-2,nverts) + 1;
        ip = mod(i,nverts) + 1;
        tm = 0; ti = 0;
        if abs(A(q,im)) > 1e-13
            tm = (sv(q,im) - D(q,im)/sv(q,i))/A(q,im);
        end
        if abs(A(q,i)) > 1e-13
            ti = (sv(q,ip) - D(q,i)/sv(q,i))/A(q,i);
        end
        w(i) = tm + ti;
        if grad_bool
            em = rrv(im,:)/sv(q,im);
            ei = rrv(i,:)/sv(q,i);
            ep = rrv(ip,:)/sv(q,ip);
            cm = em/sv(q,im) - ei/sv(q,i);
            ci = ei/sv(q,i) - ep/sv(q,ip);
            am = acos( (rrv(im,:)*rrvt(:,i)) / (sv(q,im)*sv(q,i)) );
            ai = acos( (rrv(ip,:)*rrvt(:,i)) / (sv(q,ip)*sv(q,i)) );
            tm = tan(am/2);
            ti = tan(ai/2);
            R(i,:,q) = tm/(tm+ti)*[-cm(2),cm(1)]/sin(am) + ti/(tm+ti)*[-ci(2),ci(1)]/sin(ai) + ei/sv(q,i);
        end
    end
    W = w*ont;
    bout(q,:) = w/W;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bout, R] = mv3D(bout, rv, sv, grad_bool)
[nqx, nverts] = size(sv);
R = zeros(nverts, 3, nqx);
for q=1:nqx
    if d_bool(q), continue; end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gout = get_basis_gradients(bout, R)
gout = R; gout(:,:,:) = 0;
for i=1:size(R,3)
    phi = bout(i,:);
    phiR = phi*R(:,:,i);
    for k=1:size(R,2)
        gout(:,k,i) = phi'.*(R(:,k,i) - phiR(:,k));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = get_face_R(dim, b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%