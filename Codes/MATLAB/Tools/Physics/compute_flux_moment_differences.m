%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Compute Flux Moment Differences
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, denom] = compute_flux_moment_differences(DoF, FE, flux, flux0, ngs, mom, n_type)
err = 0; denom = 0;
if isinf(n_type)
    for g=1:length(ngs)
        gg = ngs(g);
        for m=1:length(mom)
            if isempty(flux0{gg,m})
                terr = max(abs(flux{gg,m}));
                tden = max(abs(flux{gg,m}));
            else
                terr = max(abs(flux{gg,m} - flux0{gg,m}));
                tden = max(abs(flux{gg,m}));
            end
            if terr > err, err = terr; end
            if tden > denom, denom = tden; end
        end
    end
elseif isa(n_type, 'double') && n_type > 0
    % Loop through Cells and compute L2 Norm
    for c=1:DoF.TotalCells
        M = FE.CellMassMatrix{c};
        cdofs = DoF.ConnectivityArray{c};
        zn = ones(length(cdofs), 1);
        for g=1:length(ngs)
            gg = ngs(g);
            for m=1:length(mom)
                if ~isempty(flux0{gg,m})
                    err = err + (M*(flux{gg,m}(cdofs) - flux0{gg,m}(cdofs)).^n_type)'*zn;
                    denom = denom + (M*flux{gg,m}(cdofs).^n_type)'*zn;
                else
                    err = err + (M*(flux{gg,m}(cdofs)).^n_type)'*zn;
                    denom = denom + (M*flux{gg,m}(cdofs).^n_type)'*zn;
                end
            end
        end
    end
    err = sqrt(err);
    denom = sqrt(denom);
else
    error('Cannot determine norm type.');
end