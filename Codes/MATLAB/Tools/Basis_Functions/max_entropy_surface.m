%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Maximum Entropy Basis Function Generator - Surface Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the Maximum-Entropy basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = max_entropy_surface( varargin )
% Collect Input Arguments
% -----------------------
verts = varargin{1};
faces = varargin{2}; nfaces = length(faces);
flags = varargin{3};
order = varargin{4};
nverts = varargin{5};
% Prepare Vertices and Dimensional Space
% --------------------------------------
[nv,dim] = size(verts);
ns = zeros(nfaces, 1);
for f=1:nfaces
    ns(f) = length(faces{f});
end
% Quick Error Checking
% --------------------
if order > 2, error('2 is the maximum order (serendipity).'); end
if dim < 2, error('Maximum Entropy Functions only defined for 2D/3D.'); end
if sum(flags) == 0, error('No matrices specified.'); end
if sum(flags) ~= nargout, error('Insufficient output variables.'); end
% Allocate All Memory
% -------------------
[M, G] = allocate_mem(nv, nfaces, ns, dim, flags);
qx = []; qw = []; nq = zeros(nfaces, 1);
for f=1:nfaces
    [tqx, tqw] = get_general_surface_quadrature(verts, faces{f}, order + 1, false);
    qx = [qx; tqx]; qw = [qw; tqw]; nq(f) = length(tqw);
end
% Get Basis Function Values and Gradients
% ---------------------------------------
if flags(2)
    [bme, gme] = max_entropy_basis_functions(verts, qx, faces, order, nverts);
else
    bme = max_entropy_basis_functions(verts, qx, faces, order, nverts);
end
% Create Local Matrices
% ---------------------
c = 0;
for f=1:nfaces
    MM = M{f};
    fv = faces{f};
    for i=1:nq(f)
        c = c + 1;
        rqw = qw(c);
        bt = bme(c,fv);
        MM = MM + rqw * (bt'*bt);
        if flags(2)
            bt = bme(c,:);
            bg = gme(:,:,c);
            GG = G{f};
            for d=1:dim
                GG{d} = GG{d} + rqw * (bg(:,d)*bt);
            end
        end
    end
    M{f} = MM;
    if flags(2), G{f} = GG; end
end
% Set Matrix Output Arguments
% ---------------------------
counter = 1;
if flags(1) == 1
    varargout{counter} = M;
    counter = counter + 1;
end
if flags(2) == 1
    varargout{counter} = G;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, G] = allocate_mem(nv, nfaces, ns, dim, flags)
% mass
if flags(1)
    M = cell(nfaces, 1);
    for f=1:nfaces
        M{f} = zeros(ns(f), ns(f));
    end
else
    M = [];
end
% gradient
if flags(2)
    G = cell(nfaces, 1);
    znv = zeros(nv, nv);
    for f=1:nfaces
        G{f} = cell(dim, 1);
        for d=1:dim
            G{f}{d} = znv;
        end
    end
else
    G = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%