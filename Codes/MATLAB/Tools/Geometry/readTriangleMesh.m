function [nodes, elements, edges]  = readTriangleMesh(name)

elements = readele(strcat(name,'.ele'));
nodes = readnode(strcat(name,'.node'));
edges = readedge(strcat(name,'.edge'));
% neighs = readneigh(strcat(name,'.neigh'));

[m_el,n_el] = size(elements);
[m_ed,n_ed] = size(edges);
[m_n,n_n] = size(nodes);
% [m_ne,n_ne] = size(neighs);

if min(elements) == 0, elements = elements + 1; end
if min(edges(:,1:2)) == 0, edges(:,1:2) = edges(:,1:2) + 1; end
% neighs = neighs + ones(m_ne,n_ne,'int32');

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ele = readele(fname)
  fid = fopen(fname, 'rt');

  if (fid < 0)
    error('Could not open the file.');
  end;

  % read the first line of the elements file
  first_line  = fscanf(fid, '%d', 3);
  len = first_line(1);
  vert = first_line(2);
  att = first_line(3);

  ele = zeros(len, vert, 'int32');

  % Again, we use a cycle instead of matlab specific 
  %   functions like textscan, so that you can easily
  %   port to another programming language if needed.

  for i=1:len
      % read the first number on the current line
      %   (for our files this is just the element number)
      j=fscanf(fid, '%d', 1); 
%       if (j~=i-1)
%           error('Invalid file format.');
%       end

      % read the three (vert) numbers defining the
      %  vertices of the current element
      ele(i,:) = fscanf(fid, '%d', vert)';

      % finally read the attributes (which we currently discard)
      fscanf(fid, '%d', att);
  end

  fclose(fid);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function node = readnode(fname)
  fid = fopen(fname, 'rt');

  if (fid < 0)
    error('Could not open the file.');
  end;

  % read the first line of the nodes file
  first_line  = fscanf(fid, '%d', 4);
  len = first_line(1);
  dim = first_line(2);
  att = first_line(3);
  bdr = first_line(4);

  node = zeros(len, dim+bdr);

  % Again, we use a cycle instead of matlab specific 
  %   functions like textscan, so that you can easily
  %   port to another programming language if needed.

  for i=1:len
      % read the first number on the current line
      %   (for our files this is just the node number)
      j=fscanf(fid, '%d', 1); 
%       if (j~=i-1)
%           error('Invalid file format.');
%       end

      % read the two (dim) real numbers defining the
      %  coordinates of the current node
      node(i,1:dim) = fscanf(fid, '%f', dim)';

      % read the attributes (which we currently discard)
      fscanf(fid, '%d', att);

      % and the boundary marker
      if (bdr~=0)
          node(i,dim+1:dim+bdr) = fscanf(fid, '%d', bdr)';
      end
  end

  fclose(fid);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edge = readedge(fname)
  fid = fopen(fname, 'rt');

  if (fid < 0)
    error('Could not open the file.');
  end;

  % read the first line of the nodes file
  first_line  = fscanf(fid, '%d', 2);
  len = first_line(1);
  bdr = first_line(2);

  edge = zeros(len, 2+bdr);

  % Again, we use a cycle instead of matlab specific 
  %   functions like textscan, so that you can easily
  %   port to another programming language if needed.

  for i=1:len
      % read the first number on the current line
      %   (for our files this is just the edge number)
      j=fscanf(fid, '%d', 1); 
%       if (j~=i-1)
%           error('Invalid file format.');
%       end

      % read the two integers defining the end-points
      %  coordinates of the current edge
      edge(i,1:2) = fscanf(fid, '%d', 2)';

      % and the boundary marker
      if (bdr~=0)
          edge(i,3:2+bdr) = fscanf(fid, '%d', bdr)';
      end
  end

  fclose(fid);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function neigh = readneigh(fname)
  fid = fopen(fname, 'rt');

  if (fid < 0)
    error('Could not open the file.');
  end
  
  % read the first line of the neighbor file
  first_line  = fscanf(fid, '%d', 2);
  len = first_line(1);
  bdr = first_line(2);
  
  neigh = zeros(len,bdr,'int32');
  
  % Again, we use a cycle instead of matlab specific 
  %   functions like textscan, so that you can easily
  %   port to another programming language if needed.
  
  for i=1:len
      % read the first number on the current line
      %   (for our files this is just the cell number)
      j=fscanf(fid, '%d', 1);
      if (j~=i-1)
          error('Invalid file format.');
      end
      
      neigh(i,:) = fscanf(fid, '%d', bdr);
  end
  
  fclose(fid);
return