%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Orthogonal Cartesian Serendipity Main Generation Function
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB function to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the Serendipity basis functions on quad and 
%                   hex meshes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input Space:    1) Number of geometric vertices
%                   2) Vertices
%                   3) Face Vertices
%                   4) FEM Order
%                   5) Volumetric Matrix Flags
%                   6) Surface Matrix Flags
%                   7) Quadrature boolean
%                   8) Quadrature Order (Optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = bf_cell_func_Serendipity_Cartesian_Orthogonal( varargin )
% Collect Input/Output Arguments
% ------------------------------
nout = nargout;
nverts = varargin{1};
verts = varargin{2}(1:nverts,:);
faces = varargin{3}; nf = length(faces);
ord = varargin{4};
v_flags = varargin{5};
s_flags = varargin{6};
q_bool = varargin{7};
q_ord = ord+2;
if nargin > 7
    if ~isempty(varargin{8}),q_ord = varargin{8};end
end
% Prepare Vertices and Dimensional Space
% --------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
[nv,dim] = size(verts);
% Quick Error Checking
% --------------------

% Compute and exit immediately if 1D
% ------------------------------------------------------------------------------
if dim == 1
    [bf_V,bf_S,QV,QS] = bf_func_1D(varargin{:});
    varargout = {bf_V, bf_S, QV, QS};
    return
end
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxilliary Function Calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_num_dofs(dim, ord)
if dim == 2
    if ord == 1
        out = 4;
    elseif ord == 2
        out = 8;
    elseif ord == 3
        out = 12;
    end
elseif dim == 3
    if ord == 1
        out = 8;
    elseif ord == 2
        out = 20;
    elseif ord == 3
        out = 32;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_face_dofs(dim, ord)
if dim == 2
    if ord == 1
        out = [1,2;2,3;3,4;4,1];
    elseif ord == 2
        out = [1,2,5;2,3,6;3,4,7;4,1,8];
    elseif ord == 3
        out = [1,2,5,6;2,3,7,8;3,4,9,10;4,1,11,12];
    end
elseif dim == 3
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx, qw] = get_ref_quadrature(dim, ord)
[tqx, tqw] = get_legendre_gauss_quad(ord); nt = length(tqw);
if dim == 1
    qx = tqx; qw = tqw;
elseif dim == 2
    qw = tqw*tqw'; qw = qw(:); qw = qw / sum(qw);
    qxx = zeros(nt, nt); qxy = zeros(nt, nt);
    for i=1:nt
        qxx(i,:) = tqx(i);
        qxy(:,i) = tqx(i);
    end
    qx = [qxx(:), qxy(:)];
elseif dim == 3
    qw = zeros(nt, nt, nt); qqw = tqw*tqw';
    qxx = zeros(nt, nt, nt); qxy = zeros(nt, nt, nt); qxz = zeros(nt, nt, nt); 
    for i=1:nt
        qxx(i,:,:) = tqx(i);
        qxy(:,i,:) = tqx(i);
        qxz(:,:,i) = tqx(i);
        qw(:,:,i)  = qqw*tqw(i);
    end
    qx = [qxx(:), qxy(:), qxz(:)];
    qw = qw(:); qw = qw / sum(qw);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,K,G] = get_volume_matrices(dim, ord, dx, dy)
if dim == 2
    if ord==1
        M = dx*dy*[ 4, 2, 1, 2;...
                    2, 4, 2, 1;...
                    1, 2, 4, 2;...
                    2, 1, 2, 4]/36;
        K = [  (dx^2 + dy^2)/(3*dx*dy),    dx/(6*dy) - dy/(3*dx), -(dx^2 + dy^2)/(6*dx*dy),    dy/(6*dx) - dx/(3*dy);...
                 dx/(6*dy) - dy/(3*dx),  (dx^2 + dy^2)/(3*dx*dy),    dy/(6*dx) - dx/(3*dy), -(dx^2 + dy^2)/(6*dx*dy);...
              -(dx^2 + dy^2)/(6*dx*dy),    dy/(6*dx) - dx/(3*dy),  (dx^2 + dy^2)/(3*dx*dy),    dx/(6*dy) - dy/(3*dx);...
                 dy/(6*dx) - dx/(3*dy), -(dx^2 + dy^2)/(6*dx*dy),    dx/(6*dy) - dy/(3*dx),  (dx^2 + dy^2)/(3*dx*dy)];
        G{1} = dy/12 * [ -2, 2, 1, -1;...
                         -2, 2, 1, -1;...
                         -1, 1, 2, -2;...
                         -1, 1, 2, -2];
        G{2} = dx/12 * [ -2, -1, 1, 2;...
                         -1, -2, 2, 1;...
                         -1, -2, 2, 1;...
                         -2, -1, 1, 2];
    elseif ord == 2
        M = dx*dx/180 * [ 6,  2,  3,  2, -6, -8, -8, -6;...
                          2,  6,  2,  3, -6, -6, -8, -8;...
                          3,  2,  6,  2, -8, -6, -6, -8;...
                          2,  3,  2,  6, -8, -8, -6, -6;...
                         -6, -6, -8, -8, 32, 20, 16, 20;...
                         -8, -6, -6, -8, 20, 32, 20, 16;...
                         -8, -8, -6, -6, 16, 20, 32, 20;...
                         -6, -8, -8, -6, 20, 16, 20, 32];
        K = 1/90 * [ (52*(dx^2 + dy^2))/(dx*dy),    (17*dx)/dy + (28*dy)/dx, (23*(dx^2 + dy^2))/(dx*dy),    (28*dx)/dy + (17*dy)/dx,   (6*dx)/dy - (80*dy)/dx, - (40*dx)/dy - (6*dy)/dx, - (6*dx)/dy - (40*dy)/dx,   (6*dy)/dx - (80*dx)/dy;...
                        (17*dx)/dy + (28*dy)/dx, (52*(dx^2 + dy^2))/(dx*dy),    (28*dx)/dy + (17*dy)/dx, (23*(dx^2 + dy^2))/(dx*dy),   (6*dx)/dy - (80*dy)/dx,   (6*dy)/dx - (80*dx)/dy, - (6*dx)/dy - (40*dy)/dx, - (40*dx)/dy - (6*dy)/dx;...
                     (23*(dx^2 + dy^2))/(dx*dy),    (28*dx)/dy + (17*dy)/dx, (52*(dx^2 + dy^2))/(dx*dy),    (17*dx)/dy + (28*dy)/dx, - (6*dx)/dy - (40*dy)/dx,   (6*dy)/dx - (80*dx)/dy,   (6*dx)/dy - (80*dy)/dx, - (40*dx)/dy - (6*dy)/dx;...
                        (28*dx)/dy + (17*dy)/dx, (23*(dx^2 + dy^2))/(dx*dy),    (17*dx)/dy + (28*dy)/dx, (52*(dx^2 + dy^2))/(dx*dy), - (6*dx)/dy - (40*dy)/dx, - (40*dx)/dy - (6*dy)/dx,   (6*dx)/dy - (80*dy)/dx,   (6*dy)/dx - (80*dx)/dy;...
                         (6*dx)/dy - (80*dy)/dx,     (6*dx)/dy - (80*dy)/dx,   - (6*dx)/dy - (40*dy)/dx,   - (6*dx)/dy - (40*dy)/dx, (48*dx)/dy + (160*dy)/dx,                        0,  (80*dy)/dx - (48*dx)/dy,                        0;...
                       - (40*dx)/dy - (6*dy)/dx,     (6*dy)/dx - (80*dx)/dy,     (6*dy)/dx - (80*dx)/dy,   - (40*dx)/dy - (6*dy)/dx,                        0, (160*dx)/dy + (48*dy)/dx,                        0,  (80*dx)/dy - (48*dy)/dx;...
                       - (6*dx)/dy - (40*dy)/dx,   - (6*dx)/dy - (40*dy)/dx,     (6*dx)/dy - (80*dy)/dx,     (6*dx)/dy - (80*dy)/dx,  (80*dy)/dx - (48*dx)/dy,                        0, (48*dx)/dy + (160*dy)/dx,                        0;...
                         (6*dy)/dx - (80*dx)/dy,   - (40*dx)/dy - (6*dy)/dx,   - (40*dx)/dy - (6*dy)/dx,     (6*dy)/dx - (80*dx)/dy,                        0,  (80*dx)/dy - (48*dy)/dx,                        0, (160*dx)/dy + (48*dy)/dx];
        G{1} = dy/180 * [-12,  -8,  -3,   3,  20, -14,   0,  14;...
                           8,  12,  -3,   3, -20, -14,   0,  14;...
                           3,  -3,  12,   8,   0, -14, -20,  14;...
                           3,  -3,  -8, -12,   0, -14,  20,  14;...
                         -20,  20,   0,   0,   0,  40,   0, -40;...
                          14,  26,  26,  14, -40,  48, -40, -48;...
                           0,   0,  20, -20,   0,  40,   0, -40;...
                         -26, -14, -14, -26,  40,  48,  40, -48];
        
        G{2} = dx/180 * [-12,   3,  -3,  -8,  14,   0, -14,  20;...
                           3, -12,  -8,  -3,  14,  20, -14,   0;...
                           3,   8,  12,  -3,  14, -20, -14,   0;...
                           8,   3,  -3,  12,  14,   0, -14, -20;...
                         -26, -26, -14, -14, -48,  40,  48,  40;...
                           0, -20,  20,   0, -40,   0,  40,   0;...
                          14,  14,  26,  26, -48, -40,  48, -40;...
                         -20,   0,   0,  20, -40,   0,  40,   0];
    elseif ord == 3
        
    end
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%