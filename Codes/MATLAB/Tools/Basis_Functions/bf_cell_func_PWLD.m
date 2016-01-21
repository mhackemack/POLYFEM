%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          PWLD Main Generation Function
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the PWLD basis functions.
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
function varargout = bf_cell_func_PWLD( varargin )
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
    if ~isempty(varargin{8}),q_ord = varargin{9};end
end
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
% Compute and exit immediately if 1D
% ------------------------------------------------------------------------------
if dim == 1
    [bf_V,bf_S,QV,QS] = bf_cell_func_1D(varargin{:});
    varargout = {bf_V, bf_S, QV, QS};
    return
end
% Quick Error Checking
% --------------------
if ord ~= 1, error('PWLD requires 1st order FE space.'); end
% ------------------------------------------------------------------------------
% Prepare PWLD Preliminaries
% --------------------------
rcenter = mean(verts); a = 1/nv;
fcenter = get_face_centers(dim, verts, faces);
if dim == 2
    i_offset = 1;
else
    i_offset = 0;
end
[mv, ms] = get_ref_mass_matrix(dim,lump_bool);
[bv, bs, db] = get_ref_basis(dim);
nsides = get_total_sides(dim, faces);
side_areas = get_face_side_areas(dim, verts, faces);
side_vols = zeros(nsides, 1);
dbt = db';
qx_v = []; qw_v = []; bvals_v = []; bgrads_v = [];
qx_s = cell(nf,1); qw_s = cell(nf,1);
bvals_s = cell(nf,1); bgrads_s = cell(nf,1);
if q_bool
    % Volume Terms
    [rqx_v, rqw_v] = get_ref_quad(dim, q_ord);   rqw_v = rqw_v*(dim*(dim-1));
    nrqx_v = length(rqw_v); nqx_v = nsides*nrqx_v;
    qx_v = zeros(nqx_v, dim); qw_v = zeros(nqx_v, 1);
    bvals_v = zeros(nqx_v, nv); bgrads_v = zeros(nv, dim, nqx_v);
    ones_rv = ones(nrqx_v,1);
    % Surface Terms
    [rqx_s, rqw_s] = get_ref_quad(dim-1, q_ord); rqw_s = rqw_s*(dim-1);
    nrqx_s = length(rqw_s); nqx_s = nsides*nrqx_s;
    qx_s = cell(nf, 1); qw_s = cell(nf, 1);
    bvals_s = cell(nf, 1); bgrads_s = cell(nf, 1);
    ones_rs = ones(nrqx_s,1);
    % Reference Values
    if dim == 2
        b_ref = [ones(nrqx_v,1)-rqx_v(:,1)-rqx_v(:,2),rqx_v(:,1),rqx_v(:,2)];
        rqx_s = [rqx_s,zeros(nrqx_s,1)];
        b_ref_s = [ones(nrqx_s,1)-rqx_s(:,1)-rqx_s(:,2),rqx_s(:,1),rqx_s(:,2)];
    else
        b_ref = [ones(nrqx_v,1)-rqx_v(:,1)-rqx_v(:,2)-rqx_v(:,3),rqx_v(:,1),rqx_v(:,2),rqx_v(:,3)];
        rqx_s = [rqx_s,zeros(nrqx_s,1)];
        b_ref_s = [ones(nrqx_s,1)-rqx_s(:,1)-rqx_s(:,2)-rqx_s(:,3),rqx_s(:,1),rqx_s(:,2),rqx_s(:,3)];
    end
end
% Allocate Matrix Space
% ------------------------------------------------------------------------------
znv = zeros(nv);
J = zeros(dim);
M = znv;
K = znv;
G = cell(dim, 1);
for d=1:dim, G{d} = znv; end
MM = cell(nf, 1);
G2 = cell(nf, 1);
for f=1:nf
    MM{f} = zeros(length(faces{f}));
    for d=1:dim, G2{f}{d} = znv; end
end
% Collect all Matrices and Quadratures
% ------------------------------------------------------------------------------
% Loop through each face in 2D/3D
tside = 0;
for f=1:nf
    fv = faces{f}; nfv = length(fv); b = 1/nfv;
    nfe = (nfv - i_offset);
    % Loop through face edges (only 1 in 2D)
    for i=1:nfe
        tside = tside + 1;
        ii = [i,mod(i,nfv)+1];
        eee = fv(ii);
        tverts = [verts(eee,:);fcenter{f};rcenter]';
        v0 = verts(eee(1),:);
        % Calculate jacobian and additional info
        for d=1:dim
            J(:,d) = tverts(:,d+1) - tverts(:,1);
        end
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
        side_vols(tside) = svol;
        % ----------------------------------------------------------------------
        % begin 2D/3D matrix contributions
        % volume - mass matrix
        if v_flags(1)
            mm = svol*mv;
            M(eee,eee) = M(eee,eee) + mm(1:2,1:2);
            M = M + a*a*mm(end,end);
            M(eee(1),:) = M(eee(1),:) + a*mm(1,end);
            M(eee(2),:) = M(eee(2),:) + a*mm(2,end);
            M(:,eee(1)) = M(:,eee(1)) + a*mm(end,1);
            M(:,eee(2)) = M(:,eee(2)) + a*mm(end,2);
            if dim == 3
                M(fv,fv) = M(fv,fv) + b*b*mm(end-1,end-1);
                M(fv,:) = M(fv,:) + a*b*mm(end-1,end);
                M(:,fv) = M(:,fv) + a*b*mm(end,end-1);
                M(eee(1),fv) = M(eee(1),fv) + b*mm(1,end-1);
                M(eee(2),fv) = M(eee(2),fv) + b*mm(2,end-1);
                M(fv,eee(1)) = M(fv,eee(1)) + b*mm(end-1,1);
                M(fv,eee(2)) = M(fv,eee(2)) + b*mm(end-1,2);
            end
        end
        % volume - stiffness matrix
        if v_flags(2)
            s = svol*(db*(invJ*invJ')*dbt);
            K(eee,eee) = K(eee,eee) + s(1:2,1:2);
            K = K + a*a*s(end,end);
            K(eee(1),:) = K(eee(1),:) + a*s(1,end);
            K(eee(2),:) = K(eee(2),:) + a*s(2,end);
            K(:,eee(1)) = K(:,eee(1)) + a*s(end,1);
            K(:,eee(2)) = K(:,eee(2)) + a*s(end,2);
            if dim == 3
                K(fv,fv) = K(fv,fv) + b*b*s(end-1,end-1);
                K(fv,:) = K(fv,:) + a*b*s(end-1,end);
                K(:,fv) = K(:,fv) + a*b*s(end,end-1);
                K(eee(1),fv) = K(eee(1),fv) + b*s(1,end-1);
                K(eee(2),fv) = K(eee(2),fv) + b*s(2,end-1);
                K(fv,eee(1)) = K(fv,eee(1)) + b*s(end-1,1);
                K(fv,eee(2)) = K(fv,eee(2)) + b*s(end-1,2);
            end
        end
        % volume - gradient matrix
        if v_flags(3)
            c = db*invJ*detJ;
            for j=1:dim
                g = (c(:,j)*bv)';
                GG = G{j};
                GG(eee,eee) = GG(eee,eee) + g(1:2,1:2);
                GG = GG + a*a*g(end,end);
                GG(eee(1),:) = GG(eee(1),:) + a*g(1,end);
                GG(eee(2),:) = GG(eee(2),:) + a*g(2,end);
                GG(:,eee(1)) = GG(:,eee(1)) + a*g(end,1);
                GG(:,eee(2)) = GG(:,eee(2)) + a*g(end,2);
                if dim == 3
                    GG(fv,fv) = GG(fv,fv) + b*b*g(end-1,end-1);
                    GG(fv,:) = GG(fv,:) + a*b*g(end-1,end);
                    GG(:,fv) = GG(:,fv) + a*b*g(end,end-1);
                    GG(eee(1),fv) = GG(eee(1),fv) + b*g(1,end-1);
                    GG(eee(2),fv) = GG(eee(2),fv) + b*g(2,end-1);
                    GG(fv,eee(1)) = GG(fv,eee(1)) + b*g(end-1,1);
                    GG(fv,eee(2)) = GG(fv,eee(2)) + b*g(end-1,2);
                end
                G{j} = GG;
            end
        end
        % surface - mass matrix
        if s_flags(1)
            mt = side_areas(tside)*ms;
            MM{f}(ii,ii) = MM{f}(ii,ii) + mt(1:2,1:2);
            if dim == 3
                MM{f}           =  MM{f} + b*b*mt(3,3);
                MM{f}(ii(1),:) = MM{f}(ii(1),:) + b*mt(1,3);
                MM{f}(ii(2),:) = MM{f}(ii(2),:) + b*mt(2,3);
                MM{f}(:,ii(1)) = MM{f}(:,ii(1)) + b*mt(3,1);
                MM{f}(:,ii(2)) = MM{f}(:,ii(2)) + b*mt(3,2);
            end
        end
        % surface - gradient matrix
        if s_flags(2)
            c = db*invJ*side_areas(tside);
            for d=1:dim
                g = c(:,d)*bs;
                GG = G2{f}{d};
                GG(eee,eee) = GG(eee,eee) + g(1:2,1:2);
                GG(eee(1),:) = GG(eee(1),:) + a*g(1,end);
                GG(eee(2),:) = GG(eee(2),:) + a*g(2,end);
                GG(:,eee(1)) = GG(:,eee(1)) + a*g(end,1);
                GG(:,eee(2)) = GG(:,eee(2)) + a*g(end,2);
                GG = GG + a*a*g(end,end);
                if dim == 3
                    GG(fv,fv) = GG(fv,fv) + b*b*g(3,3);
                    GG(eee(1),fv)  = GG(eee(1),fv)  + b*g(1,3);
                    GG(eee(2),fv)  = GG(eee(2),fv)  + b*g(2,3);
                    GG(fv,eee(1))  = GG(fv,eee(1))  + b*g(3,1);
                    GG(fv,eee(2))  = GG(fv,eee(2))  + b*g(3,2);
                    GG(:,fv) = GG(:,fv) + a*b*g(4,3);
                    GG(fv,:) = GG(fv,:) + a*b*g(3,4);
                end
                G2{f}{d} = GG;
            end
        end
        % end 2D/3D matrix contributions
        % ----------------------------------------------------------------------
        % begin quadrature generation
        if q_bool
            tg_ref = db*invJ;
            % Volume
            iiv = (tside-1)*nrqx_v+1:tside*nrqx_v;
            qw_v(iiv) = svol*rqw_v;
            qx_v(iiv,:) = ones_rv*v0 + (J*rqx_v')';
            bvals_v(iiv,eee) = b_ref(:,1:2);
            for q=1:nrqx_v
%                 qx_v(iiv(q),:) = v0 + (J*rqx_v(q,:)')';
                % Values
%                 bvals_v(iiv(q),eee) = b_ref(q,1:2);
                bvals_v(iiv(q),:) = bvals_v(iiv(q),:) + a*b_ref(q,end);
                if dim == 3
                    bvals_v(iiv(q),fv) = bvals_v(iiv(q),fv) + b*b_ref(q,end-1);
                end
                % Gradients
                
            end
            % Surface
            tqq = ones_rs*v0 + (J*rqx_s')';
%             tqq = zeros(nrqx_s,dim); 
            tbs = zeros(nrqx_s,nfv); tgs = zeros(nv, dim, nrqx_s);
            tbs(:,ii) = b_ref_s(:,1:2);
            for q=1:nrqx_s
%                 tqq(q,:) = v0 + (J*rqx_s(q,:)')';
                % Values
%                 tbs(q,ii) = b_ref_s(q,1:2);
                if dim == 3
                    tbs(q,:) = tbs(q,:) + b*b_ref_s(q,end-1);
                end
                % Gradients
                tgs(eee,:,q) = tg_ref(1:2,:);
                for d=1:dim
%                     tgs(eee,d,q) = tg_ref(1:2,d);
                    tgs(:,d,q) = tgs(:,d,q) + a*tg_ref(end,d);
                    if dim == 3
                        tgs(fv,d,q) = tgs(fv,d,q) + b*tg_ref(end-1,d);
                    end
                end
            end
            qx_s{f} = [qx_s{f};tqq];
            qw_s{f} = [qw_s{f};side_areas(tside)*rqw_s];
            tznv = zeros(nrqx_s, nv); tznv(:,fv) = tbs;
            bvals_s{f} = [bvals_s{f};tznv];
            bgrads_s{f} = [bgrads_s{f};tgs];
        end
        % end quadrature generation
        % ----------------------------------------------------------------------
    end
end
% Process Output Structures
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G};
% Surface Matrices
varargout{2} = {MM, G2};
% Quadrature Structures
varargout{3} = {qx_v, qw_v, bvals_v, bgrads_v};
varargout{4} = {qx_s, qw_s, bvals_s, bgrads_s};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_ref_quad(dim, ord)
if dim == 1
    [qx, qw] = get_legendre_gauss_quad(ord);
elseif dim == 2
    [qx, qw] = Quad_On_Triangle(ord);
elseif dim == 3
    [qx, qw] = Quad_On_Tetra(ord);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_centers(dim, verts, faces)
nf = length(faces);
out = cell(nf, 1);
if dim == 3
    for f=1:nf
        out{f} = mean(verts(faces{f}, :));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bv, bs, db] = get_ref_basis(dim)
if dim == 2
    bv = [1, 1, 1]/6;
    bs = [1/2, 1/2, 0];
    db = [    -1    -1
               1     0
               0     1];
else
    bv = [1, 1, 1, 1]/24;
    bs = [1/3, 1/3, 1/3, 0];
    db = [    -1    -1    -1
               1     0     0
               0     1     0
               0     0     1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m, ms] = get_ref_mass_matrix(dim,lump)
if ~lump
    if dim == 2
        m = [2,1,1;1,2,1;1,1,2]./12;
        ms = [2,1,0;1,2,0;0,0,0]./6;
    else
        m = [2,1,1,1;1,2,1,1;1,1,2,1;1,1,1,2]./20;
        ms = [2,1,1,0;1,2,1,0;1,1,2,0;0,0,0,0]./12;
    end
else
    if dim == 2
        m = [1,0,0;0,1,0;0,0,1]./3;
        ms = [1,0,0;0,1,0;0,0,0]./2;
    else
        m = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1]./4;
        ms = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,0]./3;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_total_sides(dim, faces)
if dim == 2
    out = length(faces);
else
    out = 0;
    for f=1:length(faces)
        fv = faces{f};
        out = out + length(fv);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_side_areas(dim, v, faces)
if dim == 2
    out = zeros(length(faces), 1);
    for i=1:length(faces)
        out(i) = norm(diff(v(faces{i},:)));
    end
else
    out = [];
    for f=1:length(faces)
        fv = faces{f}; nfv = length(fv);
        fc = mean(v(fv,:));
        for i=1:nfv
            ii = [i, mod(i,nfv)+1];
            verts = [v(fv(ii),:);fc];
            pArea = polygonArea3d(verts);
            out = [out, pArea];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = heron_3D( fv )
a = norm(fv(2,:) - fv(1,:));
b = norm(fv(3,:) - fv(2,:));
c = norm(fv(1,:) - fv(3,:));
s = (a+b+c)/2;
out = sqrt(s*(s-a)*(s-b)*(s-c));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_grads(dim)
if dim == 2
    out = [    -1    -1
                1     0
                0     1];
else
    out = [    -1    -1    -1
                1     0     0
                0     1     0
                0     0     1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%