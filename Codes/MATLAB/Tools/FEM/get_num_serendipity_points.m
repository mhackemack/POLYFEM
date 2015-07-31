%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Serendipity Number Routine
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to return the number of  possible serendipity
%                   degrees of freedom for a 2D or 3D cell.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ntot = get_num_serendipity_points( dim, nverts, nfaces, order)
if dim == 1
    ntot = order + 1;
elseif dim == 2
    ntot = order*nverts;
elseif dim == 3
    ntot = nverts + (order-1)*(nverts+nfaces-2); % Euler's formula
end