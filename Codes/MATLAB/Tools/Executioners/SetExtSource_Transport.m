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
% Get some preliminary information
% ------------------------------------------------------------------------------
mquad = data.Quadrature(qid);
ngroups = length(groups);
XS = data.XS(xsid);
% Build the MMS source
% ------------------------------------------------------------------------------
if data.MMS.PerformMMS
    angdirs = mquad.AngularDirections;
    src = cell(ngroups,mquad.NumberAngularDirections);
    for g=1:ngroups
        for m=1:mquad.NumberAngularDirections
            src{g,m} = zeros(DoF.TotalDoFs,1);
        end
    end
    % Loop through cells
    for c=1:mesh.TotalCells
        cn = DoF.ConnectivityArray{c};
        qx = FE.CellQuadNodes{c};
        qw = FE.CellQuadWeights{c};
        cb = FE.CellBasisValues{c};
        % Loop through energy groups
        for g=1:ngroups
            gfunc = XS.ExtSource{groups(g)};
            % Loop through angular directions
            for q=1:mquad.NumberAngularDirections
                src{g,q}(cn) = src{g,q}(cn) + cb'*(qw.*gfunc(qx,angdirs(q,:)));
            end
        end
    end
end
% Build the constant distributed source
% ------------------------------------------------------------------------------
if ~data.MMS.PerformMMS
    src = cell(ngroups,1); % External sources are restricted to P0
    for g=1:ngroups
        src{g} = zeros(DoF.TotalDoFs,1);
    end
    % Loop through cells
    for c=1:mesh.TotalCells
        cmat = mesh.MatID(c);
        cn = DoF.ConnectivityArray{c};
        F = FE.CellFunctionMatrix{c};
        % Loop through energy groups
        for g=1:ngroups
            src{g,1}(cn) = src{g,1}(cn) + XS.ExtSource(cmat,groups(g))*F;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%