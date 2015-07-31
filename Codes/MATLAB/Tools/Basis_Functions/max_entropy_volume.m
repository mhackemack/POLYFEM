%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Maximum Entropy Basis Function Generator - Volume Integrals
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
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = max_entropy_volume(varargin)
% Collect Input Arguments
% -----------------------
verts = varargin{1};
faces = varargin{2};
flags = varargin{3};
order = varargin{4};
nverts = varargin{5};
% Prepare Vertices and Dimensional Space
% --------------------------------------
[mv,nv] = size(verts); 
if nv > mv, verts = verts'; end
dim = size(verts, 2);
% Quick Error Checking
% --------------------
if order > 2, error('2 is the maximum order (serendipity).'); end
if dim < 2, error('Maximum Entropy Functions only defined for 2D/3D.'); end
if sum(flags) ~= nargout, error('Insufficient output variables.'); end
% Allocate Matrix Memory
% ----------------------
ntot = get_num_serendipity_points( dim, nverts, length(faces), order);
M = zeros(ntot,ntot);       % mass matrix
K = zeros(ntot,ntot);       % stiffness matrix
G = cell(dim,1);            % gradient matrix
for i=1:dim
    G{i} = zeros(ntot,ntot);
end
% Get Basis Function Values and Gradients
% ---------------------------------------
[qx, qw] = get_general_volume_quadrature(verts, faces, 2*order+1, true); nqx = length(qw);
if flags(2) || flags(3)
    [bme, gme] = max_entropy_basis_functions(verts, qx, faces, order, nverts);
else
    bme = max_entropy_basis_functions(verts, qx, faces, order, nverts);
end
% Create Local Matrices
% ---------------------
% mass matrix
for q=1:nqx
    bt = bme(q,:);
    M = M + qw(q) * (bt'*bt);
end
% stiffness matrix
if flags(2)
    for q=1:nqx
        bg = gme(:,:,q);
        K = K + qw(q) * (bg*bg');
    end
end
% gradient matrix
if flags(3)
    for q=1:nqx
        bt = bme(q,:);
        bg = gme(:,:,q);
        for d=1:dim
            G{d} = G{d} + qw(q) * (bg(:,d)*bt)';
        end
    end
end
% Set Matrix Output Arguments
% ---------------------------
counter = 1;
if flags(1)
    varargout{counter} = M;
    counter = counter + 1;
end
if flags(2)
    varargout{counter} = K;
    counter = counter + 1;
end
if flags(3)
    varargout{counter} = G;
end

