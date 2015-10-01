%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Set Boundary Incident Fluxes
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
function data = set_incident_boundary_fluxes(data,xsid,qid,mesh,DoF)
global glob
% Set some initial information and define memory space
% ------------------------------------------------------------------------------
mquad = data.Quadrature(qid);
nangs = mquad.NumberAngularDirections;
XS = data.XS(xsid);
ngroups = data.Groups.NumberEnergyGroups;
data.Fluxes.IncomingBoundaryFlux = cell(mesh.TotalFaces,1);
% Loop through boundary faces and set incident fluxes
% ------------------------------------------------------------------------------
for ff=1:mesh.TotalBoundaryFaces
    % Set cell structure for face
    f = mesh.BoundaryFaces(ff);
    data.Fluxes.IncomingBoundaryFlux{f} = cell(nangs,ngroups);
    % Retrieve some face information
    fnorm = mesh.FaceNormal(f,:)';
    fid = mesh.FaceID(f);
    % Allocate beam flux if necessary - then skips to next face afterwards
    if XS.BCFlags(fid) == glob.IncidentBeam
        % Allocate some info
        min_dot = 0.0; min_ind = 1;
        nfcn = length(DoF.FaceCellNodes{f,1});
        z = ones(nfcn,1);
        % Determine Maximum Angle
        for m=1:nangs
            adir = mquad.AngularDirections(m,:);
            fdot = adir*fnorm;
            if fdot < 0 && fdot < min_dot && abs(fdot - min_dot) > 1e-12
                min_dot = fdot; min_ind = m;
            end
        end
        % Loop through energies and set fluxes
        for g=1:ngroups
            val = XS.BCVals{fid}(g);
            data.Fluxes.IncomingBoundaryFlux{f}{min_ind,g} = val*z;
        end
        % Skip to next face
        continue
    end
    % Loop through angular directions and find incoming directions
    for m=1:nangs
        adir = mquad.AngularDirections(m,:);
        fdot = adir*fnorm;
        if fdot < 0
            fcn = DoF.FaceCellNodes{f,1}; nfcn = length(fcn);
            % Set incident fluxes by boundary condition type
            if XS.BCFlags(fid) == glob.Vacuum || XS.BCFlags(fid) == glob.Reflecting
                z = zeros(nfcn,1);
                for g=1:ngroups
                    data.Fluxes.IncomingBoundaryFlux{f}{m,g} = z;
                end
            elseif XS.BCFlags(fid) == glob.IncidentIsotropic
                z = ones(nfcn,1);
                bcv = XS.BCVals{fid};
                for g=1:ngroups
                    data.Fluxes.IncomingBoundaryFlux{f}{m,g} = bcv(g)*z/mquad.AngQuadNorm;
                end
            elseif XS.BCFlags(fid) == glob.IncidentCurrent
                % never actually did this one...
            elseif XS.BCFlags(fid) == glob.Function
                bcv = XS.BCVals{fid};
                fxn = DoF.NodeLocations(fcn,:);
                for g=1:ngroups
                    fvals = bcv{g}(fxn,adir);
                    data.Fluxes.IncomingBoundaryFlux{f}{m,g} = fvals;
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%