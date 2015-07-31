%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Calculate Reflecting Angle
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = calculate_reflecting_angles( ndat, mesh )
global glob
dim = mesh.Dimension;
% Calculate 1D Reflecting Angles
if dim == 1
    ndat = calculate_1D_reflecting_angles(ndat, mesh);
    return
end
na = ndat.Transport.NumberAngularDirections;
octAngles = na / 2^dim;
octNumbers = get_octant_numbers( na, mesh.Dimension );
ndat.Transport.ReflectingBoundaryAngles = cell(mesh.TotalFaces, 1);
mins = get_mins(dim);
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fid = mesh.FaceID(ff);
    nflag = ndat.Transport.BCFlags(fid);
    if nflag ~= glob.Reflecting, continue; end
    ndat.Transport.ReflectingBoundaryAngles{ff} = zeros(na, 1);
    fnorm = mesh.FaceNormal(ff,:);
    for f_out=1:na
        angDir = ndat.Transport.AngularDirections(f_out,:);
        fdot = angDir * fnorm';
        if abs(fdot) < glob.small, continue; end
        if fdot > 0
            f_in = get_reflected_angle(fnorm, f_out, octNumbers(f_out), octAngles, mins);
            ndat.Transport.ReflectingBoundaryAngles{ff}(f_out) = f_in;
            ndat.Transport.ReflectingBoundaryAngles{ff}(f_in) = f_out;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = calculate_1D_reflecting_angles(ndat, mesh)
global glob
ndat.Transport.ReflectingBoundaryAngles = cell(mesh.TotalFaces, 1);
na = ndat.Transport.NumberAngularDirections;
% Left Boundary
fid = mesh.FaceID(1);
nflag = ndat.Transport.BCFlags(fid);
if nflag == glob.Reflecting
    ndat.Transport.ReflectingBoundaryAngles{1} = zeros(na, 1);
    c = na;
    for m=1:na/2
        ndat.Transport.ReflectingBoundaryAngles{1}(m) = c;
        ndat.Transport.ReflectingBoundaryAngles{1}(c) = m;
        c = c - 1;
    end
end
% Right Boundary
fid = mesh.FaceID(end);
nflag = ndat.Transport.BCFlags(fid);
if nflag == glob.Reflecting
    ndat.Transport.ReflectingBoundaryAngles{end} = zeros(na, 1);
    c = na/2;
    for m=na/2+1:na
        ndat.Transport.ReflectingBoundaryAngles{end}(m) = c;
        ndat.Transport.ReflectingBoundaryAngles{end}(c) = m;
        c = c - 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_mins(dim)
if dim == 1
    out = get_1D_mins();
elseif dim == 2
    out = get_2D_mins();
else
    out = get_3D_mins();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_1D_mins()
out.xmin = -1;
out.xmax =  1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_2D_mins()
out.xmin = [-1,  0];
out.xmax = [ 1,  0];
out.ymin = [ 0, -1];
out.ymax = [ 0,  1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_3D_mins()
out.xmin = [-1,  0,  0];
out.xmax = [ 1,  0,  0];
out.ymin = [ 0, -1,  0];
out.ymax = [ 0,  1,  0];
out.zmin = [ 0,  0, -1];
out.zmax = [ 0,  0,  1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_octant_numbers( na, dim )
out = zeros( na, 1 );
octAngles = na / 2^dim;
for m=1:na
    if m <= octAngles
        out(m) = 1;
    elseif m > octAngles && m <= 2*octAngles
        out(m) = 2;
    elseif m > 2*octAngles && m <= 3*octAngles
        out(m) = 3;
    elseif m > 3*octAngles && m <= 4*octAngles
        out(m) = 4;
    elseif m > 4*octAngles && m <= 5*octAngles
        out(m) = 5;
    elseif m > 5*octAngles && m <= 6*octAngles
        out(m) = 6;
    elseif m > 6*octAngles && m <= 7*octAngles
        out(m) = 7;
    elseif m > 7*octAngles && m <= 8*octAngles
        out(m) = 8;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_reflected_angle( fnorm, m , octant, octAngles, mins )
global glob
% xmin
if norm(fnorm - mins.xmin) < glob.small
    if octant == 2
        out = m - octAngles;
    elseif octant == 3
        out = m + octAngles;
    elseif octant == 6
        out = m - octAngles;
    elseif octant == 7
        out = m + octAngles;
    else
        out = m;
    end
    return
end
% xmax
if norm(fnorm - mins.xmax) < glob.small
    if octant == 1
        out = m + octAngles;
    elseif octant == 4
        out = m - octAngles;
    elseif octant == 5
        out = m + octAngles;
    elseif octant == 8
        out = m - octAngles;
    else
        out = m;
    end
    return
end
% ymin
if norm(fnorm - mins.ymin) < glob.small
    if octant == 3
        out = m - octAngles;
    elseif octant == 4
        out = m - 3*octAngles;
    elseif octant == 7
        out = m - octAngles;
    elseif octant == 8
        out = m - 3*octAngles;
    else
        out = m;
    end
    return
end
% ymax
if norm(fnorm - mins.ymax) < glob.small
    if octant == 1
        out = m + 3*octAngles;
    elseif octant == 2
        out = m + octAngles;
    elseif octant == 5
        out = m + 3*octAngles;
    elseif octant == 6
        out = m + octAngles;
    else
        out = m;
    end
    return
end
% zmin
if norm(fnorm - mins.zmin) < glob.small
    if octant == 5
        out = m - 4*octAngles;
    elseif octant == 6
        out = m - 4*octAngles;
    elseif octant == 7
        out = m - 4*octAngles;
    elseif octant == 8
        out = m - 4*octAngles;
    else
        out = m;
    end
    return
end
% zmax
if norm(fnorm - mins.zmax) < glob.small
    if octant == 1
        out = m + 4*octAngles;
    elseif octant == 2
        out = m + 4*octAngles;
    elseif octant == 3
        out = m + 4*octAngles;
    elseif octant == 4
        out = m + 4*octAngles;
    else
        out = m;
    end
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%