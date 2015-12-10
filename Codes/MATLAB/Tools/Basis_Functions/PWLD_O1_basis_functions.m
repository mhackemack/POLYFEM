%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to produce the basis function values and
%                   gradients for the linear case.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PWLD_O1_basis_functions(varargin)
% Collect Input Arguments
% ------------------------------------------------------------------------------
verts  = varargin{1};
qx     = varargin{2};
faces  = varargin{3};
order  = varargin{4};
nverts = varargin{5};
grad_bool = false;
% Determine Input Characteristics and Perform Error Checking
% ------------------------------------------------------------------------------
[nv, dim] = size(verts); nfaces = length(faces);
nqx = size(qx, 1);
if nargout > 1, grad_bool = true; end
if order > 1, error('Only 1st order in this functor. Go fix your code.'); end
if nv ~= nverts, error('Number of vertices does not align. Go fix your code.'); end
% Allocate Memory Space
% ------------------------------------------------------------------------------
nvones = ones(1,nv);
bout = zeros(nqx, nv); gout = zeros(nv, dim, nqx);
% Retrieve 1D Values and Gradients
% ------------------------------------------------------------------------------
if dim == 1
    x0 = verts(1);
    J = get_simplex_jacobian(verts,dim); invJ = 1/J;
    xref = invJ*(qx - x0);
    bout = get_ref_values(dim, xref);
    if grad_bool
        refgrad = get_ref_grads(dim); trg = refgrad * invJ;
        for q=1:nqx
            gout(:,:,q) = trg;
        end
    end
end
% Retrieve 2D Values and Gradients
% ------------------------------------------------------------------------------
if dim == 2
    rcenter = mean(verts); q_nums = (1:nqx)';
    if grad_bool, refgrad = get_ref_grads(dim);  end
    % Loop through faces
    for f=1:nfaces
        % Build edge triangle
        fv = faces{f}; fverts = verts(fv,:);
        vv = [fverts;rcenter];
        % Build triangle jacobian
        v0 = vv(1,:); J = get_simplex_jacobian(vv, dim);
        detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
        invJ = [J(2,2),-J(1,2);-J(2,1),J(1,1)]/detJ;
        % Determine points inside triangle
        qxbool = inpoly(qx, vv);
        tq_nums = q_nums(qxbool);
        xtemp = qx(tq_nums,:); nxt = size(xtemp,1);
        % Map to the reference triangle and get reference values
        v0ref = ones(nxt,1)*v0; xref = (invJ*(xtemp - v0ref)')';
        tbvals = get_ref_values(dim, xref);
        bout(tq_nums,fv) = bout(tq_nums,fv) + tbvals(:,1:2);
        bout(tq_nums,:) = bout(tq_nums,:) + (1/nv)*tbvals(:,end)*nvones;
        if grad_bool
            tqg = refgrad * invJ; 
            tvqg = (1/nv)*nvones'*tqg(end,:);
            for q=1:nxt
                gout(fv,:,tq_nums(q)) = gout(fv,:,tq_nums(q)) + tqg(1:2,:);
                gout(:,:,tq_nums(q)) = gout(:,:,tq_nums(q)) + tvqg;
            end
        end
    end
    bout = bout./(sum(bout,2)*nvones);
end
% Retrieve 3D Values and Gradients
% ------------------------------------------------------------------------------
if dim == 3
    rcenter = mean(verts); q_nums = (1:nqx)';
    if grad_bool, refgrad = get_ref_grads(dim);  end
    % Loop through faces
    for f=1:nfaces
        fv = faces{f}; nf = length(fv);
        fverts = verts(fv,:); fcenter = mean(fverts);
        % Loop through face vertices
        for ff=1:nf
            fi = [ff, mod(ff,nf)+1];
        end
    end
end
% Assign Output Arguments
% ------------------------------------------------------------------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = get_simplex_jacobian(verts, dim)
J = zeros(dim);
vverts = verts';
for d=1:dim
    J(:,d) = vverts(:,d+1) - vverts(:,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_values(dim, xin)
if dim == 1
    out = [1-xin,xin];
elseif dim == 2
    x = xin(:,1); y = xin(:,2);
    out = [1-x-y, x, y];
else
    x = xin(:,1); y = xin(:,2); z = xin(:,3);
    out = [1-x-y-z, x, y, z];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_grads(dim)
if dim == 1
    out = [-1;1];
elseif dim == 2
    out = [-1,-1;1,0;0,1];
else
    out = [-1,-1,-1;1,0,0;0,1,0;0,0,1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%