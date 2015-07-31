%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Interpolate Refinement Solution
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
%   Notes:          As of right now, the interpolation is handled by the simple
%                   average of flux values from 1 cell refinement level up using
%                   the average geometric distances between nodes. This acts to
%                   get the solution 'close' to the refined solution for the
%                   next level so that we do not have to recompute the entire
%                   solution space.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d2,flux] = interpolate_ref_flux(d1, d2, mesh, dof1, dof2, flux0, flux)
ng = d2.Neutronics.numberEnergyGroups;
nm = d2.Neutronics.TotalFluxMoments;
% Loop through refined mesh cells
for c=1:mesh.TotalCells
    c0 = mesh.PreviousCell(c);
    cdofs1 = dof1.ConnectivityArray{c0}; ncdofs1 = length(cdofs1);
    cdofs2 = dof2.ConnectivityArray{c};  ncdofs2 = length(cdofs2);
    cn1 = dof1.NodeLocations(cdofs1,:); cn2 = dof2.NodeLocations(cdofs2,:);
    wts = zeros(ncdofs2, ncdofs1);
    % Loop through old/new dofs
    for i=1:ncdofs2
        for j=1:ncdofs1
            len = norm(cn2(i,:)-cn1(j,:));
            if len < 1e-12
                wts(i,:) = 0;
                wts(i,j) = 1;
                break;
            end
            wts(i,j) = 1/len;
        end
    end
    sumwts = sum(wts,2);
    % Loop through DoFs
    for i=1:ncdofs2
        % Loop through all energy groups
        for g=1:ng
            % Loop through all flux moments
            for m=1:nm
                f0 = flux0{g,m}(cdofs1);
                flux{g,m}(cdofs2) = (wts*f0)./sumwts;
            end
        end
    end
end