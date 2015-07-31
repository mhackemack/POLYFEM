%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Wachpress Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the basis function and gradient
%                   values using the Wachpress Methodology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = wachspress_basis_functions( varargin )
nout = nargout;
grad_bool = false;
% Collect Input Arguments
% -----------------------
verts = varargin{1};
qx = varargin{2};
faces = varargin{3};
% Prepare Vertices and Dimensional Space
% --------------------------------------
dim = size(verts, 2);
% Quick Error Checking
% --------------------
% if order > 1, error('Wachpress only defined for order 1.'); end
% Allocate Matrix Memory
% ----------------------
if nout > 1, grad_bool = true; end
if dim == 2
    [b, bg] = wach2d( verts, qx, grad_bool );
elseif dim == 3
    [b, bg] = wach3d( verts, faces, qx, grad_bool );
end
% Set Output Arguments
% --------------------
varargout{1} = b;
if grad_bool, varargout{2} = bg; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b, bg] = wach2d( v, qx, grad_bool )
global glob
n = size(qx, 1); nv = size(v, 1);
w = zeros(n, nv); R = zeros(nv, 2, n); p = zeros(nv, 2, n);
un = get_2D_normals( v ); unt = un';
b = zeros(n, nv);
face_bool = logical(zeros(n,1));
for i=1:n
    x = qx(i,:);
    for j=1:nv
        h = (v(j,:) - x)*unt(:,j);
        if abs(h) < glob.small
            face_bool(i) = true;
        else
            p(j,:,i) = un(j,:) / h;
        end
    end
end
for i=1:n
    if ~face_bool(i), continue; end
    x = qx(i,:);
    for j=1:nv
        jj = mod(j,nv) + 1;
        sj = (v(j,:) - x); hj = (v(j,:) - x)*unt(:,j);
        sjj = (v(jj,:) - x); hjj = (v(jj,:) - x)*unt(:,j);
        if abs(hj) < glob.small && abs(hjj) < glob.small
            sj = norm(sj); sjj = norm(sjj);
            rr = sj + sjj;
            b(i,j) = sjj/rr;
            b(i,jj) = sj/rr;
        end
    end
end
for i=1:n
    for j=1:nv
        jm1 = mod(j-2,nv) + 1;
        w(i,j) = det([p(jm1,:,i);p(j,:,i)]);
        R(j,:,i) = p(jm1,:,i) + p(j,:,i);
    end
    if ~face_bool(i), b(i,:) = w(i,:) / sum(w(i,:)); end
end
% Get gradients
bg = get_gradient_terms( b, R, n, nv, 2, grad_bool );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b, bg] = wach3d( v, f, qx, grad_bool )
n = size(qx, 1); nv = size(v, 1); nf = length(f);
w = zeros(n, nv); R = zeros(nv, 3, n); p = zeros(nv, 3, n);
b = zeros(n, nv);
for i=1:n
    
    for j=1:nv
        
    end
    b(i,:) = w(i,:) / sum(w(i,:));
end
% Get gradients
bg = get_gradient_terms( b, R, n, nv, 3, grad_bool );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_2D_normals( v )
n = size(v,1);
out = zeros(n,2);
for i=1:n
    d = v(mod(i,n)+1,:) - v(i,:);
    out(i,:) = [d(2),-d(1)]/norm(d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_gradient_terms( b, R, n, nv, dim, grad_bool )
if ~grad_bool
    out = []; return
end
out = zeros(nv, dim, n);
for i=1:n
    phi = b(i,:);
    phiR = phi*R(:,:,i);
    for k=1:dim
        out(:,k,i) = phi'.*(R(:,k,i) - phiR(:,k));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_neighbor_graphs( v, f )
[nv, dim] = size(v);
nf = length(f);
out = cell(nv, 1);

for i=1:nv
    vt = [];
    for j=1:nf
        nfv = length(f{j});
        for k=1:nfv
            
        end
    end
    % form unique set
    out({i}) = unique(vt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%