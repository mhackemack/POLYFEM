%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process AMR Data
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = process_amr_data(data)
if ~isfield(data,'AMR')
    data.AMR.RefineMesh = false;
elseif ~isfield(data.AMR, 'refineMesh')
    data.AMR.RefineMesh = false;
elseif data.AMR.RefineMesh
    if ~isfield(data.AMR,'refinementLevels')
        error('# of refinement levels required.')
    elseif data.AMR.RefinementLevels < 1
        data.AMR.refineMesh = false;
        return
    end
    if ~isfield(data.AMR, 'AMRIrregularity')
        data.AMR.AMRIrregularity = inf;
    else
        if ~isnumeric(data.AMR.AMRIrregularity)
            error('Cannot determine AMR regularity');
        else
            data.AMR.AMRIrregularity = round(data.AMR.AMRIrregularity);
        end
    end
    if ~isfield(data.AMR, 'ProjectSolution')
        data.AMR.ProjectSolution = 0;
    end
    if data.AMR.RefinementTolerance < 0 && abs(data.AMR.RefinementTolerance) > 1e-13
        error('Refinement tolerance needs to be between 0 and 1.');
    else
        data.AMR.RefinementTolerance = 0.0;
    end
    if data.AMR.RefinementTolerance > 1 && abs(data.AMR.RefinementTolerance - 1) > 1e-13
        error('Refinement tolerance needs to be between 0 and 1.');
    else
        data.AMR.RefinementTolerance = 1.0;
    end
    if ~isfield(data.AMR, 'refinementType')
        error('Need to specify a refinement type.')
    end
end