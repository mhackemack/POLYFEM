%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Determine thermal groups
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
function out = determine_thermal_groups(S)
% Restrict to 0th moments
S = S(:,:,1); ng = size(S,1);
% Loop through energy groups and find thermal cutoff
gcutoff = 0;
for g=1:ng
    if abs(sum(S(g,g+1:end))) > 1e-12
        gcutoff = g;
        break
    end
end
% Write thermal groups to list
out = gcutoff:ng;