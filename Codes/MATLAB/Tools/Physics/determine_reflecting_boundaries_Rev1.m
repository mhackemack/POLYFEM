%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Determine Reflecting Boundaries (Rev1)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = determine_reflecting_boundaries_Rev1( data, mesh )
global glob
dim = mesh.Dimension;
xsid = data.Transport.XSID;
mins = get_geometry_mins(mesh);
data.Transport.HasReflectingBoundary = false;
data.Transport.HasOpposingReflectingBoundary = false;
data.Transport.OpposingReflectingBoundaryDimension = false(mesh.Dimension, 1);
data.Transport.ReflectingBoundaries = false(mesh.Dimension, 2);
for ff=1:mesh.TotalBoundaryFaces
    f = mesh.BoundaryFaces(ff);
    fid = mesh.FaceID(f);
    fc = mesh.FaceCenter(f,:);
    if data.XS(xsid).BCFlags(fid) == glob.Reflecting
        data.Transport.HasReflectingBoundary = true;
        for d=1:dim
            % Check min
            if abs(fc(d) - mins(d,1)) < 1e-13
                data.Transport.ReflectingBoundaries(d,1) = true;
            end
            % Check max
            if abs(fc(d) - mins(d,2)) < 1e-13
                data.Transport.ReflectingBoundaries(d,2) = true;
            end
        end
    end
end
% Determine Opposing Reflecting Boundaries
for d=1:dim
    if data.Transport.ReflectingBoundaries(d,1) && data.Transport.ReflectingBoundaries(d,2)
        data.Transport.HasOpposingReflectingBoundary = true;
        data.Transport.OpposingReflectingBoundaryDimension(d) = true;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_geometry_mins(mesh)
dim = mesh.Dimension;
out = [mesh.minX, mesh.maxX];
if dim > 1, out = [out;[mesh.minY, mesh.maxY]]; end
if dim > 2, out = [out;[mesh.minZ, mesh.maxZ]]; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%