%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Check Angular Quadrature
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
function data = process_acceleration_data(data)
% Set some Acceleration information if field is empty
if ~isfield(data,'Acceleration') || ~data.Transport.PerformAcceleration
    ngs = data.Groups.NumberGroupSets;
    data.Acceleration.WGSAccelerationBool = false(ngs,1);
    data.Acceleration.AGSAccelerationBool = false;
    data.Acceleration.WGSAccelerationResidual = false(ngs,1);
    data.Acceleration.AGSAccelerationResidual = false;
    data.Acceleration.WGSAccelerationID = zeros(ngs,1);
    data.Acceleration.AGSAccelerationID = 0;
    % Set Eigenvalue acceleration values
    if strcmpi(data.problem.ProblemType, 'eigenvalue')
        data.Acceleration.PIAccelerationBool = false;
    end
    return
end
% Perform some additional error checking
