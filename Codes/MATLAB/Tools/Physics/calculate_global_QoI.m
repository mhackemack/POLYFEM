%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate Global QoI
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to calculate the global QoI.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = calculate_global_QoI(data, xsid, mesh, DoF, FE, flux, varargin)
% Process Input Information
% ------------------------------------------------------------------------------
nm = data.problem.NumberMaterials;
if iscell(flux)
    ng = size(flux,1);
else
    ng = 1;
end
if nargin > 5
    qt = get_qoi_type( data.XS(xsid), varargin{1}, nm, ng );
else
    qt = get_qoi_type( data.XS(xsid), 'Flux', nm, ng );
end
% Loop through mesh cells
out = 0;
for c=1:mesh.TotalCells
    matID = mesh.MatID(c);
    dc = DoF.ConnectivityArray{c};
    M = FE.CellMassMatrix{c};
    for g=1:ng
        out = out + qt(matID,g)*sum(M*flux{g,1}(dc));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_qoi_type( XS, s, nm, ng )
switch(s)
    case('Total')
        out = XS.TotalXS;
    case('Absorption')
        out = XS.AbsorbXS;
    case('Fission')
        out = XS.FissionXS;
    case('Production')
        out = XS.FissionXS.*XS.NuBar;
    case('Flux')
        out = ones(nm, ng);
    otherwise
        out = ones(nm, ng);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%