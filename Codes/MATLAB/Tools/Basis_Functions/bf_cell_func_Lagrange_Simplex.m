%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Simplex Lagrange Main Generation Function
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB function to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the Lagrange basis functions on simplexes 
%                   (triangles and tetrahedra).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input Space:    1) Number of geometric vertices
%                   2) Vertices
%                   3) Face Vertices
%                   4) FEM Order
%                   5) FEM Lumping Boolean
%                   6) Volumetric Matrix Flags
%                   7) Surface Matrix Flags
%                   8) Quadrature boolean
%                   9) Quadrature Order (Optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = bf_cell_func_Lagrange_Simplex( varargin )
% Collect Input/Output Arguments
% ------------------------------------------------------------------------------
nout = nargout;
nverts = varargin{1};
verts = varargin{2}(1:nverts,:);
faces = varargin{3}; nf = length(faces);
ord = varargin{4};
lump_bool = varargin{5};
v_flags = varargin{6};
s_flags = varargin{7};
q_bool = varargin{8};
q_ord = ord+2;
if nargin > 8
    if ~isempty(varargin{9}),q_ord = varargin{9};end
end
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
% Compute and exit immediately if 1D
% ------------------------------------------------------------------------------
if dim == 1
    [bf_V,bf_S,QV,QS] = bf_func_1D(varargin{:});
    varargout = {bf_V, bf_S, QV, QS};
    return
end
% Quick Error Checking
% ------------------------------------------------------------------------------
if dim == 2
    if nv > 4 || nv < 3, error('Only triangles and quads in 2D.'); end
    if ord > 3, error('Cannot go higher than 3rd order for 2D at this time.'); end
end
if dim == 3
    if ord > 1, error('Cannot go higher than 1st order for 3D at this time.'); end
end
% Generate Basis Function Space
% ------------------------------------------------------------------------------
% Retrieve Reference Function Handles
% -----------------------------------
ref_val_func = get_lagrange_function('vals', 1, dim);
ref_grad_func = get_lagrange_function('grads', 1, dim);
% Generate Reference Quadrature Spaces
% ------------------------------------
if dim == 2
    ref_x = [0,0;1,0;0,1];
    [rqx_V, rqw_V] = Quad_On_Triangle(q_ord); rqw_V = 2*rqw_V;
    [rqx_S, rqw_S] = get_legendre_gauss_quad(q_ord);
elseif dim == 3
    ref_x = [0,0,0;1,0,0;0,1,0;0,0,1];
    [rqx_V, rqw_V] = Quad_On_Tetra(q_ord); rqw_V = 6*rqw_V;
    [rqx_S, rqw_S] = Quad_On_Triangle(q_ord); rqw_S = 2*rqw_S;
end
nq_V = length(rqw_V);
nq_S = length(rqw_S);
zs = ones(nq_S, 1);
% Compute Jacobian
% ----------------
J = zeros(dim);
x0 = verts(1,:); tverts = verts';
for d=1:dim
    J(:,d) = tverts(:,d+1) - tverts(:,1);
end
% We directly compute the Jacobian determinants and inverses to drastically save
% on code execution times for larger problems.
if dim==2
    detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
    invJ = [J(2,2),-J(1,2);-J(2,1),J(1,1)]/detJ;
    svol = detJ/2;
else
    detJ = J(1,1)*(J(3,3)*J(2,2)-J(3,2)*J(2,3)) - J(2,1)*(J(3,3)*J(1,2)-J(3,2)*J(1,3)) + J(3,1)*(J(2,3)*J(1,2)-J(2,2)*J(1,3));
    invJ = [(J(3,3)*J(2,2)-J(3,2)*J(2,3)),-(J(3,3)*J(1,2)-J(3,2)*J(1,3)), (J(2,3)*J(1,2)-J(2,2)*J(1,3));...
           -(J(3,3)*J(2,1)-J(3,1)*J(2,3)), (J(3,3)*J(1,1)-J(3,1)*J(1,3)),-(J(2,3)*J(1,1)-J(2,1)*J(1,3));...
            (J(3,2)*J(2,1)-J(3,1)*J(2,2)),-(J(3,2)*J(1,1)-J(3,1)*J(1,2)), (J(2,2)*J(1,1)-J(2,1)*J(1,2))]/detJ;
    svol = detJ/6;
end
% Compute Real Space Quadrature and Basis Values/Gradients
% --------------------------------------------------------
ntot = get_num_dofs(dim, ord); zt = ones(1,ntot);
f_dofs = get_face_dofs(dim, ord);
% Volume
qx_V = zeros(nq_V, dim); qw_V = rqw_V*svol;
bvals_v = ref_val_func(ord, rqx_V);
bgrads_v = zeros(ntot,dim,nq_V);
for q=1:nq_V
    qx_V(q,:) = x0 + (J*rqx_V(q,:)')';
    tg = ref_grad_func(ord, rqx_V(q,:));
    bgrads_v(:,:,q) = tg*invJ;
end
% Surface
bvals_s = cell(nf, 1); bgrads_s = cell(nf, 1);
qx_S    = cell(nf, 1); qw_S     = cell(nf, 1);
for f=1:nf
    bgrads_s{f} = zeros(ntot, dim, nq_S);
    fv = faces{f};
    fverts = verts(fv,:);
    f0 = fverts(1,:);
    if dim == 2
        dx = diff(fverts);
        len = norm(dx);
        qw_S{f} = rqw_S*len;
        qx_S{f} = zs*f0 + rqx_S*dx;
        % Get ref quad then basis values/gradients
        ff = [f,mod(f,nf)+1];
        ff0 = ref_x(ff(1),:);
        dff = diff(ref_x(ff,:));
        tq = zs*ff0 + rqx_S*dff;
        bvals_s{f} = ref_val_func(ord,tq);
        for q=1:nq_S
            bgrads_s{f}(:,:,q) = ref_grad_func(ord, tq(q,:))*invJ;
        end
    elseif dim == 3
        
    end
end
% Build Matrices
% --------------
M = bvals_v'*(bvals_v.*(qw_V*zt));
if lump_bool, M = diag(sum(M)); end
K  = zeros(ntot);
G  = cell(dim,1);
for d=1:dim
    G{d} = zeros(ntot);
end
MM = cell(nf,1);
G2 = cell(nf,1);
F  = cell(nf, 1);
% Stiffness Matrix
if v_flags(2)
    for q=1:nq_V
        bg = bgrads_v(:,:,q);
        K = K + qw_V(q) * (bg*bg');
    end
end
% Gradient Matrix
if v_flags(3)
    for q=1:nq_V
        bt = bvals_v(q,:);
        bg = bgrads_v(:,:,q);
        for d=1:dim
            G{d} = G{d} + qw_V(q) * (bg(:,d)*bt)';
        end
    end
end
% Surface Matrices
for f=1:nf
    ffdd = f_dofs(f,:);
    bv = bvals_s{f};
    tmf = bv'*(bv.*(qw_S{f}*zt));
    MM{f} = tmf(ffdd,ffdd);
    if s_flags(2)
        G2{f} = cell(dim,1);
        for d=1:dim
            bg = squeeze(bgrads_s{f}(:,d,:));
            G2{f}{d} = bg*(bv.*(qw_S{f}*zt));
        end
    end
end
% Assign Output Arguments
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G};
% Surface Matrices
varargout{2} = {MM, G2, F};
% Quadrature Structure - Volume
varargout{3} = {qx_V, qw_V, bvals_v, bgrads_v};
% Quadrature Structure - Surface
varargout{4} = {qx_S, qw_S, bvals_s, bgrads_s};
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxilliary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_num_dofs(dim, ord)
if dim == 2
    if ord == 1
        out = 3;
    elseif ord == 2
        out = 6;
    elseif ord == 3
        out = 10;
    end
elseif dim == 3
    if ord == 1
        out = 4;
    elseif ord == 2
        out = 10;
    elseif ord == 3
        out = 21; % maybe
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_dofs(dim, ord)
if dim == 2
    if ord == 1
        out = [1,2;2,3;3,1];
    elseif ord == 2
        out = [1,2,4;2,3,5;3,1,6];
    elseif ord == 3
        out = [1,2,4,5;2,3,6,7;3,1,8,9];
    end
elseif dim == 3
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%