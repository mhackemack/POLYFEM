%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Average Material QoI Solution
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to calculate the average QoI over all materials.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = calculate_average_material_QoI( data, mesh, DoF, FE, flux, varargin)
% Process Input Information
% -------------------------
nm = data.problem.NumberMaterials;
if iscell(flux)
    [ng, nf] = size(flux);
else
    ng = 1; nf = 1;
end
if nargin > 5
    qt = get_qoi_type( data, varargin{1}, nm, ng );
else
    qt = get_qoi_type( data, 'Flux', nm, ng );
end
% Allocate memory
qoi = zeros(nm, ng, nf);
% Loop through mesh cells
totvol = zeros(nm, 1);
for c=1:mesh.TotalCells
    % Cell spatial info
    matID = mesh.MatID(c);
    cv = mesh.CellVolume(c);
    totvol(matID) = totvol(matID) + cv;
    dc = DoF.ConnectivityArray{c};
    M = FE.CellMassMatrix{c};
    for g=1:ng
        for f=1:nf
            qoi(matID,g,f) = qoi(matID,g,f) + qt(matID,g)*sum(M*flux{g,f}(dc));
        end
    end
end
% Divide through by material volume
for m=1:nm
    qoi(m,:,:) = qoi(m,:,:)./totvol(m);
end
% Set outputs
varargout{1} = qoi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_qoi_type( data, s, nm, ng )
if strcmp(data.Neutronics.transportMethod,'Transport')
    nd = data.Neutronics.Transport;
elseif strcmp(data.Neutronics.transportMethod,'Diffusion')
    nd = data.Neutronics.Diffusion;
end
switch(s)
    case('Total')
        out = nd.TotalXS;
    case('Absorption')
        out = nd.AbsorbXS;
    case('Fission')
        out = nd.FissionXS;
    case('Production')
        out = nd.FissionXS.*nd.NuBar;
    case('Flux')
        out = ones(nm, ng);
    otherwise
        out = ones(nm, ng);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%