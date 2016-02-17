%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Build Function Handle
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
function b_func = get_build_function(data)
TM = data.Neutronics.TransportMethod;
TT = data.Neutronics.Transport.transportType;
DM = data.Neutronics.DSAType;
if ~data.Neutronics.PerformAcceleration
    b_func = str2func(['func_build_',TM,'_',TT]);
else
    b_func = str2func(['func_build_',TM,'_',TT,'_',DM]);
end