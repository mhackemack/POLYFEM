%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Discontinuous (PWLD) Basis Function 
%                   Generator - Surface
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the global basis function 
%                   values and gradients.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PWLD_basis_functions_surface(varargin)
% Collect Input Arguments
% -----------------------
verts = varargin{1};
xin   = varargin{2};
fnodes = varargin{3};
% Determine Input Characteristics
% -------------------------------
[nv, dim] = size(verts);
nqx = size(xin, 1);
if nargout == 2
    grad_bool = true;
    ref_grads = get_ref_grads(dim);
else
    grad_bool = false;
end
if dim == 2
    bool_2D = true;
elseif dim == 3
    bool_2D = false;
end
% Get Cell Information
% --------------------
nf = length(fnodes);
rcenter = mean(verts);
a = 1 / nv;
dim_denom = dim*(dim-1);
% Allocate Matrices
% -----------------
q_vals  = zeros(nqx, nv);
q_grads = zeros(nv,dim,nqx);
q_nums  = (1:nqx)';
% Switch between dimension
% ------------------------
fverts = verts(fnodes,:);
if dim == 2
    vv = [fverts;rcenter];
    v0 = vv(1,:);
    xref = zeros(nqx,dim);
    J = get_simplex_jacobian(vv, dim);
    invJ = inv(J);
    for q=1:nqx
        xref(q,:) = invJ*(xin(q,:) - v0)';
    end
    % values
    tbvals = get_ref_values(dim, xref);
    q_vals(:,fnodes) = q_vals(:,fnodes) + tbvals(:,1:2);
    q_vals(:,:) = q_vals(:,:) + a*tbvals(:,end)*ones(1,nv);
    % gradients
    if grad_bool
        for q=1:nqx
            tqg = ref_grads * invJ;
            q_grads(fnodes,:,q) = q_grads(fnodes,:,q) + tqg(1:2,:);
            q_grads(:,:,q) = q_grads(:,:,q) + a*ones(nv,1)*tqg(end,:);
        end
    end
else
    tet_order = get_tet_ordering();
    fcenter = mean(fverts); 
    bb = 1/nf; 
    for j=1:nf
        if j==nf
            jj = [j,1];
        else
            jj = [j,j+1];
        end
        ee = fnodes(jj);
        vv = [fverts(jj,:);fcenter;rcenter];
        xbool = inpolyhedron(tet_order, vv, qx);
        v0 = vv(1,:); tq_nums = q_nums(xbool);
        xtemp = xin; xtemp(~xbool,:) = []; nxt = size(xtemp,1); xref = zeros(nxt,dim);
        J = get_simplex_jacobian(vv, dim);
        invJ = inv(J);
        for q=1:nxt
            xref(q,:) = invJ*(xtemp(q,:) - v0)';
        end
        % values
        tbvals = get_ref_values(dim, xref)/dim_denom;
        q_vals(tq_nums,ee) = q_vals(tq_nums,ee) + tbvals(:,1:2);
        q_vals(tq_nums,:) = q_vals(tq_nums,:) + a*tbvals(:,end)*ones(1,nv);
        q_vals(tq_nums,ff) = q_vals(tq_nums,ff) + bb*tbvals(:,3)*ones(1,nf);
        % gradients
        if grad_bool
            for q=1:nxt
                tqg = ref_grads * invJ;
                q_grads(ff(jj),:,tq_nums(q)) = q_grads(ff(jj),:,tq_nums(q)) + tqg(1:2,:);
                q_grads(:,:,tq_nums(q)) = q_grads(:,:,tq_nums(q)) + a*ones(nv,1)*tqg(end,:);
                q_grads(ff,:,tq_nums(q)) = q_grads(ff,:,tq_nums(q)) + bb*ones(nf,1)*tqg(3,:);
            end
        end
    end
end
% Set Outputs
% -----------
for q=1:nqx
    q_vals(q,:) = q_vals(q,:) / sum(q_vals(q,:));
end
varargout{1} = q_vals;
if grad_bool
    varargout{2} = q_grads;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = get_simplex_jacobian(verts, dim)
J = zeros(dim);
vverts = verts';
for d=1:dim
    J(:,d) = vverts(:,d+1) - vverts(:,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_values(dim, xin)
if dim == 2
    x = xin(:,1); y = xin(:,2);
    out = [1-x-y, x, y];
else
    x = xin(:,1); y = xin(:,2); z = xin(:,3);
    out = [1-x-y-z, x, y, z];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_grads(dim)
if dim == 2
    out = [-1,-1;1,0;0,1];
else
    out = [-1,-1,-1;1,0,0;0,1,0;0,0,1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_tet_ordering()
out = [1,2,3;1,2,4;1,3,4;2,3,4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%