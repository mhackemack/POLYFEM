function [nodes, eles, faces, edges]  = readTetgenMesh(name)

eles = read_tet_file(name,'.ele');
nodes = read_tet_file(name,'.node');
faces = read_tet_file(name,'.face');
edges = read_tet_file(name,'.edge');

% Increment indices by 1 if zero-indexing was used in mesh generation
if min(min(eles)) == 0
    eles = eles + 1;
end
if min(min(faces(:,1:end-1))) == 0
    faces(:,1:end-1) = faces(:,1:end-1) + 1;
end
if min(min(edges(:,1:end-1))) == 0
    edges(:,1:end-1) = edges(:,1:end-1) + 1;
end

% convert file types to integers
eles = uint32(eles);
faces = uint32(faces);
edges = uint32(edges);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = read_tet_file(fname,ename)
% get output info and allocate memory
[len,col,att,fline,desig] = get_file_info(fname,ename);
out = zeros(len,col);
% open the file
fullname = strcat(fname,ename);
fid = fopen(fullname, 'rt');
discard = fscanf(fid, '%d', fline);
% loop through file and get attributes
for i=1:len
    j=fscanf(fid, '%d', 1);
    out(i,1:end-att) = fscanf(fid, desig, col-att)';
    
    if logical(att)
        j = fscanf(fid, '%d', att);
        out(i,end) = j;
    end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [len,col,att,fline,desig] = get_file_info(fname,ename)
if strcmp(ename, '.ele')
    fline = 3;
    desig = '%d';
elseif strcmp(ename, '.node')
    fline = 4;
    desig = '%f';
elseif strcmp(ename, '.face') || strcmp(ename, '.edge')
    fline = 2;
    desig = '%d';
end

fullname = strcat(fname,ename);
fid = fopen(fullname, 'rt');
if (fid < 0), error('Could not open the file.'); end

f_line  = fscanf(fid, '%d', fline);
len = f_line(1);

if strcmp(ename, '.ele')
    att = f_line(3);   col = 4 + att;
elseif strcmp(ename, '.node')
    att = f_line(3);   col = 3 + att;
elseif strcmp(ename, '.face')
    att = f_line(2);   col = 3 + att;
elseif strcmp(ename, '.edge')
    att = f_line(2);   col = 2 + att;
end
fclose(fid);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
