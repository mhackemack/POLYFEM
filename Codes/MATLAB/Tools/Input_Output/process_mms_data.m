%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Process MMS Input Data
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
function data = process_mms_data(data)
if ~isfield(data,'MMS')
    data.MMS.PerformMMS = 0;
elseif ~isfield(data.MMS,'PerformMMS')
    data.MMS.PerformMMS = 0;
elseif data.MMS.PerformMMS
    % Set a minimum quadrature order if necessary
    if ~isfield(data.MMS,'QuadOrder')
        data.MMS.QuadOrder = 6;
        warning('Setting MMS quadrature order to 6.');
    elseif isnumeric(data.MMS.QuadOrder)
        if data.MMS.QuadOrder < 1
            data.MMS.QuadOrder = 6;
            warning('Setting MMS quadrature order to 6.');
        end
    else
        data.MMS.QuadOrder = 6;
        warning('Setting MMS quadrature order to 6.');
    end
    % Check Exact Solutions
    for g=1:data.Groups.NumberEnergyGroups
        if ~isa(data.MMS.ExactSolution{g},'function_handle')
            error('Solution entry for energy group %d is not a function handle.',g);
        end
    end
end