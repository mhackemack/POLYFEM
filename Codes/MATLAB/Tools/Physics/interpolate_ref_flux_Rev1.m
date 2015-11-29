%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Interpolate Refinement Solution (Rev1)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to interpolate a flux solution from a
%                   coarser mesh onto a finer mesh only.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:          The interpolation is performed 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flux = interpolate_ref_flux_Rev1(mesh, dof1, dof2, fe1, flux0)
[ng, nm] = size(flux0);
% Build outgoing flux structure
flux = cell(ng, nm); ndof = dof2.TotalDoFs;
for g=1:ng
    for m=1:nm
        flux(g,m) = zeros(ndof,1);
    end
end
% Loop through refined mesh cells
for c=1:mesh.TotalCells
    % Continue if mesh cell was not refined
    if ~mesh.CellRefinedLastCycle(c)
        cdofs1 = dof1.ConnectivityArray{c};
        cdofs2 = dof2.ConnectivityArray{c};
        % Loop through energy groups and moments
        for g=1:ng
            for m=1:nm
                flux{g,m}(cdofs2) = flux0{g,m}(cdofs1);
            end
        end
        continue
    end
    % Get degree of freedom information between refinements
    c0 = mesh.PreviousCell(c);
    cdofs1 = dof1.ConnectivityArray{c0}; ncdofs1 = length(cdofs1);
    cdofs2 = dof2.ConnectivityArray{c};  ncdofs2 = length(cdofs2);
    
end