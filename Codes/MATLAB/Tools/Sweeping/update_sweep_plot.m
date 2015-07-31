%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Update Sweep Plot
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_sweep_plot(mesh, cellID)
dim = mesh.Dimension;
if dim == 1
    
elseif dim == 2
    cv = mesh.CellVerts{cellID};
    v = mesh.Vertices(cv,:);
    patch(v(:,1), v(:,2), 'b');
elseif dim == 3
    
end
pause(.05)