%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          0th-order Scattering Matrices
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
function S0 = build_0th_scattering_matrices(lam, input)
% Copy Input Space
% ------------------------------------------------------------------------------
data = input.data;
mesh = input.mesh;
dof = input.dof;
fe = input.fe;
% Retrieve Preliminary Data
% ------------------------------------------------------------------------------
ndofs = dof.TotalDoFs;

% Loop through Cells and Build Matrices
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    
end