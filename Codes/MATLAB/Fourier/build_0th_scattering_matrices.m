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
function S0 = build_0th_scattering_matrices(data,mesh,dof,fe)
% Retrieve Preliminary Data
% ------------------------------------------------------------------------------
ng = data.Neutronics.numberEnergyGroups; ndofs = dof.TotalDoFs;
sxs = data.Neutronics.Transport.ScatteringXS;
S0 = zeros(ndofs,ndofs,ng,ng);
% Loop through Cells and Build Matrices
% ------------------------------------------------------------------------------
for c=1:mesh.TotalCells
    cn  = dof.ConnectivityArray{c};
    mat = mesh.MatID(c);
    M   = fe.CellMassMatrix{c};
    for g=1:ng
        for gg=1:ng
            S0(cn,cn,g,gg) = sxs(mat,g,gg)*M;
        end
    end
end