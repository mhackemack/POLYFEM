%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set External Transport Source
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function src = SetExtSource_Transport(data,xsid,qid,groups,mesh,DoF,FE)



% Build the MMS source
% ------------------------------------------------------------------------------
if data.MMS.PerformMMS
    % Loop through cells
    for c=1:mesh.TotalCells
        
    end
end
% Build the constant distributed source
% ------------------------------------------------------------------------------
if ~data.MMS.PerformMMS
    % Loop through cells
    for c=1:mesh.TotalCells
        
    end
end