%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Retrieve Cell Vertex Numbers
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
function out = retrieve_cell_vertex_numbers(sol)
n = length(sol); nmax = 0;
if n == 1
    out = sol.CellVertexNumbers;
    return
end
% Determine Maximum Vertex Number
for i=1:n
    tnm = length(sol{i}.CellVertexNumbers);
    if tnm > nmax, nmax = tnm; end
end
% Retrieve Numbers
out = zeros(n,nmax);
for i=1:n
    tnm = length(sol{i}.CellVertexNumbers);
    out(i,1:tnm) = sol{i}.CellVertexNumbers;
end