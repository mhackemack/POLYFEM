%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          2D Lagrange Basis Function Generator - Volume Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the Triangular and Quadrilateral basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = Lagrange_surface(verts, fverts, flags, deg, nverts)
nin = nargin;
nout = nargout;
%
% Collect Input Arguments
% -----------------------

%
% Perform Error Checking Operations
% ---------------------------------
[m,n] = size(verts);
if n>m, verts = verts'; end
[nv,dim] = size(verts);
if dim == 1
    if nverts > 2, error('Only 2 points in 1D.'); end
    ctype = 1;
elseif dim == 2
    if nverts ~= 3 && nverts ~= 4, error(''); end
    if nverts == 3, ctype = 1; end
    if nverts == 4, ctype = 2; end
elseif dim == 3
    if nverts ~= 4 && nverts ~= 8, error(''); end
    if nverts == 4, ctype = 1; end
    if nverts == 8, ctype = 2; end
end
%
% Get Jacobian and Other Information
% ----------------------------------
nbf = get_num_vol_basis_funcs(ctype, dim, deg);
if dim == 1
    sArea = 1.0;
elseif dim == 2
    sArea = norm(diff(verts(fverts,:)));
    [qqx, qw] = get_legendre_gauss_quad(deg+1);
	[qx, ndim] = get_2D_surface_ref_terms(ctype, fverts, qqx);
    qqx = qx;
    nqx = length(qw);
elseif dim == 3
    sArea = polygonArea3d(verts(fverts,:));
    [qx, qqx, qw, ndim] = get_3D_surface_ref_terms(verts, fverts, ctype, deg+1);
    nqx = length(qw);
end
if ctype == 1
    [~, invJ, ~, ~] = get_jacobian(ctype, verts, dim);
    if flags(2)
        db = evaluate_ref_gradients(ctype, dim, deg, qx, nqx, nbf);
        b = evaluate_ref_vals(ctype, dim, deg, qx);
    end
else
    J = evaluate_numerical_jacobian(dim, verts, fverts, qqx, qx, ndim, deg);
    invJ = cell(nqx,1);
    if dim == 3
        detJ = zeros(nqx,1);
        for q=1:nqx
            invJ{q} = inv(J{q});
            detJ(q) = det(J{q});
        end
    else
        detJ = sArea*ones(nqx,1);
        for q=1:nqx
            invJ{q} = inv(J{q});
        end
    end
    db = evaluate_ref_gradients(ctype, dim, deg, qx, nqx, nbf);
    b = evaluate_ref_vals(ctype, dim, deg, qx);
end
% 
% Switch Cell Type and Degree
% ---------------------------
counter = 1;
% mass matrix
if flags(1)
    if ctype == 1
        varargout{counter} = sArea*get_ref_mass_matrix(ctype, dim, deg);
    else
        M = zeros(nbf, nbf);
        for q=1:nqx
            M = M + (b(q,:)'*b(q,:))*qw(q)*detJ(q);
        end
%         sbf = get_surf_dofs(ctype, dim, deg, fverts);
%         varargout{counter} = M(sbf,sbf);
        varargout{counter} = M(fverts,fverts);
    end
    counter = counter + 1;
end
% gradient matrix
if flags(2)
    G = cell(dim, 1);
    for d=1:dim 
        G{d} = zeros(nbf, nbf);
    end
    switch(ctype)
        case(1)
            tJ = invJ*sArea;
            for q=1:nqx
                tdb = db(:,:,q)*tJ;
                for d=1:dim
                    G{d} = G{d} + (tdb(:,d)*b(q,:))*qw(q);
                end
            end
        case(2)
            for q=1:nqx
                tdb = db(:,:,q)*invJ{q}*detJ(q);
                for d=1:dim
                    G{d} = G{d} + (tdb(:,d)*b(q,:))*qw(q);
                end
            end
    end
    varargout{counter} = G;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, invJ, detJ, vol] = get_jacobian(ctype, verts, dim)
if ctype == 1
    J = zeros(dim);
    vverts = verts';
    for d=1:dim
        J(:,d) = vverts(:,d+1) - vverts(:,1);
    end
    invJ = inv(J);
    detJ = det(J);
    if dim == 1
        vol = detJ;
    else
        vol = detJ / (dim * (dim - 1));
    end
elseif ctype == 2
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = get_ref_mass_matrix(ctype, dim, deg)
if ctype == 1
    if dim == 2
        if deg == 1
            m = [2,1;1,2]./6;
        elseif deg == 2
            m = [4,-1,2;-1,4,2;2,2,16]./30;
        end
    elseif dim == 3
        
    end
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xout, ndim] = get_2D_surface_ref_terms(ctype, fverts, pts_1D)
if ctype == 1
    switch(fverts(1))
        case(1)
            x0 = [0,0]; dx = [1,0]; dl = 1.0;
        case(2)
            x0 = [1,0]; dx = [-1,1] / norm([-1,1]);  dl = sqrt(2)/2;
            %qw = qw/sqrt(2);
        case(3)
            x0 = [0,1]; dx = [0,-1]; dl = 1.0;
    end
    ndim = 0;
else
    switch(fverts(1))
        case(1)
            x0 = [0,0]; dx = [1,0];  dl = 1.0; ndim = 2;
        case(2)
            x0 = [1,0]; dx = [0,1];  dl = 1.0; ndim = 1;
        case(3)
            x0 = [1,1]; dx = [-1,0]; dl = 1.0; ndim = 2;
        case(4)
            x0 = [0,1]; dx = [0,-1]; dl = 1.0; ndim = 1;
    end
end
nx = length(pts_1D);
nones = ones(nx,1);
xout = [x0(1)*nones+dx(1)*pts_1D/dl, x0(2)*nones+dx(2)*pts_1D/dl];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_numerical_jacobian(dim, verts, fverts, x, xx, ind, deg)
n = size(xx,1);
out = cell(n,1);
% Evaluate Jacobian
if dim == 2
%     db = eval_quad_grads_for_jac(deg, x);
%     for i=1:n
%         out{i} = [db(fverts,1,i)'*verts(fverts,1), db(fverts,2,i)'*verts(fverts,1);...
%                   db(fverts,1,i)'*verts(fverts,2), db(fverts,2,i)'*verts(fverts,2)];
%     end
    
    xall = verts(1,1) - verts(2,1) + verts(3,1) - verts(4,1);
    yall = verts(1,2) - verts(2,2) + verts(3,2) - verts(4,2);
    x21 = verts(2,1) - verts(1,1); y21 = verts(2,2) - verts(1,2);
    x41 = verts(4,1) - verts(1,1); y41 = verts(4,2) - verts(1,2);
    for i=1:n
        out{i} = [x21 + xall*x(i,2), x41 + xall*x(i,1);y21 + yall*x(i,2), y41 + yall*x(i,1)];
    end
elseif dim == 3
%     db = eval_quad_grads_for_jac(x);
%     vverts = verts(fverts,:);
%     for i=1:n
%         J = [db(:,1,i)'*vverts(:,1), db(:,2,i)'*vverts(:,1);...
%              db(:,1,i)'*vverts(:,2), db(:,2,i)'*vverts(:,2);...
%              db(:,1,i)'*vverts(:,3), db(:,2,i)'*vverts(:,3)];
%         tc = cross([db(:,1,i)'*vverts(:,1),db(:,1,i)'*vverts(:,2),db(:,1,i)'*vverts(:,3)], [db(:,2,i)'*vverts(:,1),db(:,2,i)'*vverts(:,2),db(:,2,i)'*vverts(:,3)]);
%         out{i} = J*J';
%     end
    
    db = eval_hex_grads_for_jac(deg, xx);
    for i=1:n
        out{i} = [db(fverts,1,i)'*verts(fverts,1), db(fverts,2,i)'*verts(fverts,1), db(fverts,3,i)'*verts(fverts,1);...
                  db(fverts,1,i)'*verts(fverts,2), db(fverts,2,i)'*verts(fverts,2), db(fverts,3,i)'*verts(fverts,2);...
                  db(fverts,1,i)'*verts(fverts,3), db(fverts,2,i)'*verts(fverts,3), db(fverts,3,i)'*verts(fverts,3)];
        
        
%         out{i} = [db(:,1,i)'*verts(:,1), db(:,2,i)'*verts(:,1), db(:,3,i)'*verts(:,1);...
%                   db(:,1,i)'*verts(:,2), db(:,2,i)'*verts(:,2), db(:,3,i)'*verts(:,2);...
%                   db(:,1,i)'*verts(:,3), db(:,2,i)'*verts(:,3), db(:,3,i)'*verts(:,3)];
    end
%     if abs(sum(verts(fverts,ind))) < 1e-14
        for i=1:n
            out{i}(ind,ind) = 1.0;
        end
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_num_vol_basis_funcs(ctype, dim, deg)
if dim == 1
    out = 2 + (deg-1);
elseif dim == 2
    if ctype == 1
        if deg == 1 || deg == 2
            out = 3*deg;
        elseif deg == 3
            out = 10;
        end
    else
        out = (deg+1)^2;
    end
else
    if ctype == 1
        if deg == 1
            out = 4;
        elseif deg == 2
            out = 10;
        end
    else
        out = (deg+1)^3;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_num_surf_basis_funcs(ctype, dim, deg)
if dim == 1
    out = 1;
elseif dim == 2
    out = 2 + (deg-1);
else
    if ctype == 1
        if deg == 1 || deg == 2
            out = 3*deg;
        elseif deg == 3
            out = 10;
        end
    else
        out = (deg+1)^2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_surf_dofs(ctype, dim, deg, fverts)
if dim == 1
    out = fverts;
    return
elseif dim == 2
    f0 = fverts(1);
    v0 = 2+ctype;
    out = [fverts, v0+(deg-1)*(f0-1)+1:v0+(deg-1)*(f0-1)+(deg-1)];
else
    v0 = 4*ctype;
    if ctype == 1
        
    else
        if deg == 1
            out = fverts;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_ref_vals(ctype, dim, deg, qx)
fun = get_basis_funcs(ctype, dim);
out = fun(deg, qx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_ref_gradients(ctype, dim, deg, qx, nqx, n)
out = zeros(n,dim,nqx);
fun = get_basis_grad_funcs(ctype, dim);
for q=1:nqx
    out(:,:,q) = fun(deg, qx(q,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_funcs(ctype, dim)
out = get_lagrange_function('vals', ctype, dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_grad_funcs(ctype, dim)
out = get_lagrange_function('grads', ctype, dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qqx, qw, ndim] = get_3D_surface_ref_terms(verts, fverts, ctype, deg)
if ctype == 1
    [qqx, qw] = Quad_On_Triangle(6);
else
    [qqx, qqw] = get_legendre_gauss_quad(deg+1); nt = length(qqw);
    qw = qqw*qqw'; qw = qw(:); qw = qw / sum(qw);
    qxx = zeros(nt, nt); qxy = zeros(nt, nt);
    for i=1:nt
        qxx(i,:) = qqx(i);
        qxy(:,i) = qqx(i);
    end
    qqx = [qxx(:), qxy(:)];
    qx = zeros(length(qw), 3);
    [cdim, ndim, dval] = get_hex_face_info(fverts);
    qx(:,cdim) = qqx; qx(:,ndim) = dval;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = eval_hex_grads_for_jac(deg, xx)
n = size(xx,1); out = zeros((deg+1)^3,3,n);
for i=1:n
x = xx(i,1); y = xx(i,2); z = xx(i,3);
out(:,:,i) = [  -y.*z+y+z-1,  -x.*z+x+z-1,  -x.*y+x+y-1;...
                (1-y).*(1-z),  x.*z-x,       x.*y-x;...
                y.*(1-z),      x.*(1-z),    -x.*y;...
                y.*z-y,       (1-x).*(1-z),  x.*y-y;...
                y.*z-z,       x.*z-z         (1-x).*(1-y);...
                (1-y).*z,     -x.*z          x.*(1-y);...
                y.*z,          x.*z,         x.*y;...
                -y.*z,         (1-x).*z,     (1-x).*y];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = eval_quad_grads_for_jac(deg, xx)
n = size(xx,1); out = zeros((deg+1)^2,2,n);
for i=1:n
    s = xx(i,1); t = xx(i,2);
    if deg == 1
        out(:,:,i) = [t-1,s-1;1-t,-s;t,s;-t,1-s];
    elseif deg == 2
        out(:,:,i) = [  8*s*t*t-12*s*t+4*s-6*t*t+9*t-3,    8*s*s*t-6*s*s-12*s*t+9*s+4*t-3;...
                 8*s*t*t-12*s*t+4*s-2*t*t+3*t-1,    8*s*s*t-6*s*s-4*s*t+3*s;...
                 8*s*t*t-4*s*t-2*t*t+t,             8*s*s*t-2*s*s-4*s*t+s;...
                 8*s*t*t-4*s*t-6*t*t+3*t,           8*s*s*t-2*s*s-12*s*t+3*s+4*t-1;...
               -16*s*t*t+24*s*t-8*s+8*t*t-12*t+4, -16*s*s*t+12*s*s+16*s*t-12*s;...
               -16*s*t*t+16*s*t+4*t*t-4*t,        -16*s*s*t+8*s*s+8*s*t-4*s;...
               -16*s*t*t+8*s*t+8*t*t-4*t,         -16*s*s*t+4*s*s+16*s*t-4*s;...
               -16*s*t*t+16*s*t+12*t*t-12*t,      -16*s*s*t+8*s*s+24*s*t-12*s-8*t+4;...
                32*s*t*t-32*s*t-16*t*t+16*t,       32*s*s*t-16*s*s-32*s*t+16*s];
    elseif deg == 3

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cdim, ndim, dval] = get_hex_face_info(fverts)
tf = sort(fverts);
if isequal(tf,[1,2,3,4])
    cdim = [1,2]; ndim = 3; dval = 0;
elseif isequal(tf,[5,6,7,8])
    cdim = [1,2]; ndim = 3; dval = 1.0;
elseif isequal(tf,[1,2,5,6])
    cdim = [1,3]; ndim = 2; dval = 0.0;
elseif isequal(tf,[2,3,6,7])
    cdim = [2,3]; ndim = 1; dval = 1.0;
elseif isequal(tf,[3,4,7,8])
    cdim = [1,3]; ndim = 2; dval = 1.0;
elseif isequal(tf,[1,4,5,8])
    cdim = [2,3]; ndim = 1; dval = 0.0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%