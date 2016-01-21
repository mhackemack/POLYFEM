%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the elementary volume and
%                   surface matrices, along with the appropriate quadrature
%                   set outputs for the 1D basis functions.
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
function varargout = bf_cell_func_1D( varargin )
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
% Quick Error Checking
% ------------------------------------------------------------------------------
if nverts > 2, error('Only 2 cell vertices for 1D input.'); end
% Prepare Vertices and Dimensional Space
% ------------------------------------------------------------------------------
if size(verts,2) > 1, verts = verts'; end
[nv, dim] = size(verts);
nf = length(faces);
ntot = 2 + (ord-1);
h = verts(2) - verts(1);
xx = [verts', verts(1) + (1:(ord-1))*h/(ord)];
% Construct Volume Matrix Space
% ------------------------------------------------------------------------------
[qx, qw] = lgwt(q_ord, verts(1), verts(2)); 
zt = ones(1,ntot);
% Retrieve Basis Function Values/Gradients
bvals_v  = get_1D_values(ord, xx, qx); 
bgrads_v = get_1D_gradients(ord, xx, qx);
% Build Volumetric Matrices
M = bvals_v'*(bvals_v.*(qw*zt));
K = bgrads_v'*(bgrads_v.*(qw*zt));
G{1} = (bgrads_v'*(bvals_v.*(qw*zt)))';

% Construct Surface Matrix Space
% ------------------------------------------------------------------------------
MM = cell(nf,1);
G2 = cell(nf,1);
bvals_s  = get_1D_values(ord, xx, verts); 
bgrads_s = get_1D_gradients(ord, xx, verts);
for f=1:nf
    MM{f} = 1.0;
    G2{f}{1} = bgrads_s(f,:)'*bvals_s(f,:);
end

% Process Output Structures
% ------------------------------------------------------------------------------
% Volume Matrices
varargout{1} = {M, K, G};
% Surface Matrices
varargout{2} = {MM, G2};
% Quadrature Structure - Volume
varargout{3} = {qx, qw, bvals_v, bgrads_v};
% Quadrature Structure - Surface
varargout{4} = {{verts(1),verts(2)},{1,1},{bvals_s(1,:),bvals_s(2,:)},{bgrads_s(1,:)',bgrads_s(2,:)'}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_1D_values(ord, v, x)
if ord == 1
    out = [(v(2)-x), (x-v(1))]/(v(2) - v(1));
elseif ord == 2
    out = [(x-v(2)).*(x-v(3))/(v(1)-v(2))/(v(1)-v(3)),...
           (x-v(1)).*(x-v(3))/(v(2)-v(1))/(v(2)-v(3)),...
           (x-v(1)).*(x-v(2))/(v(3)-v(1))/(v(3)-v(2))];
elseif ord == 3
    out = [(x-v(2)).*(x-v(3)).*(x-v(4))/(v(1)-v(2))/(v(1)-v(3))/(v(1)-v(4)),...
           (x-v(1)).*(x-v(3)).*(x-v(4))/(v(2)-v(1))/(v(2)-v(3))/(v(2)-v(4)),...
           (x-v(1)).*(x-v(2)).*(x-v(4))/(v(3)-v(1))/(v(3)-v(2))/(v(3)-v(4)),...
           (x-v(1)).*(x-v(2)).*(x-v(3))/(v(4)-v(1))/(v(4)-v(2))/(v(4)-v(3))];
else
    out = ones(length(x),ord+1);
    for i=1:ord+1
        for j=1:ord+1
            if i==j, continue; end
            out(:,i) = out(:,i).*(x-v(j))/(v(i)-v(j));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_1D_gradients(ord, v, x)
if ord == 1
    out = ones(length(x),1) * [-1, 1]/(v(2) - v(1));
elseif ord == 2
    out = [  (x - v(2))/((v(1) - v(2))*(v(1) - v(3))) + (x - v(3))/((v(1) - v(2))*(v(1) - v(3))), ...
           - (x - v(1))/((v(1) - v(2))*(v(2) - v(3))) - (x - v(3))/((v(1) - v(2))*(v(2) - v(3))), ...
             (x - v(1))/((v(1) - v(3))*(v(2) - v(3))) + (x - v(2))/((v(1) - v(3))*(v(2) - v(3)))];
elseif ord == 3
    v1 = v(1); v2 = v(2); v3 = v(3); v4 = v(4);
    out = [   ((v2 - x).*(v3 - x))/((v1 - v2)*(v1 - v3)*(v1 - v4)) + ((v2 - x).*(v4 - x))/((v1 - v2)*(v1 - v3)*(v1 - v4)) + ((v3 - x).*(v4 - x))/((v1 - v2)*(v1 - v3)*(v1 - v4)),...
            - ((v1 - x).*(v3 - x))/((v1 - v2)*(v2 - v3)*(v2 - v4)) - ((v1 - x).*(v4 - x))/((v1 - v2)*(v2 - v3)*(v2 - v4)) - ((v3 - x).*(v4 - x))/((v1 - v2)*(v2 - v3)*(v2 - v4)),...
              ((v1 - x).*(v2 - x))/((v1 - v3)*(v2 - v3)*(v3 - v4)) + ((v1 - x).*(v4 - x))/((v1 - v3)*(v2 - v3)*(v3 - v4)) + ((v2 - x).*(v4 - x))/((v1 - v3)*(v2 - v3)*(v3 - v4)),...
            - ((v1 - x).*(v2 - x))/((v1 - v4)*(v2 - v4)*(v3 - v4)) - ((v1 - x).*(v3 - x))/((v1 - v4)*(v2 - v4)*(v3 - v4)) - ((v2 - x).*(v3 - x))/((v1 - v4)*(v2 - v4)*(v3 - v4))];
else
    nx = length(x); klist = (1:ord+1)';
    out = zeros(nx,ord+1); tx = ones(nx,1);
    for i=1:ord+1
        tlist = klist; tlist(i) = []; denom = 1;
        for k=1:ord
            ttx = tx;
            for kk=1:ord
                if k==kk, continue; end
                ttx = ttx.*(x-v(tlist(kk)));
            end
            out(:,i) = out(:,i) + ttx;
        end
        % Calculate denominator
        for j=1:ord
            jj = tlist(j);
            denom = denom*(v(i) - v(jj));
        end
        out(:,i) = out(:,i) / denom;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%