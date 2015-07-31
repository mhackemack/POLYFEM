function facedemo(n)

% FACETEST: Example polygonal geometries for MESHFACES.
%
%  facedemo(n)
%
% Darren Engwirda - 2007

switch n
   case 1

      node = [0.0, 0.0; 1.0,0.0; 1.0,1.0; 0.0,1.0; 1.01,0.0; 1.01,1.0; 3.0,0.0; 3.0,1.0];
      edge = [1,2; 2,3; 3,4; 4,1; 2,5; 5,6; 6,3; 5,7; 7,8; 8,6];
      face{1} = [1,2,3,4];
      face{2} = [5,6,7,2];
      face{3} = [8,9,10,6];

      meshfaces(node,edge,face);

   case 2

      % Geometry
      dtheta = pi/3;
      theta = (-pi:dtheta:(pi-dtheta))';
      node1 = [1.2*cos(theta), 1.2*sin(theta)];
      node2 = [-2.0,-2.0; 2.0,-2.0; 2.0,2.0; -2.0, 2.0];
      edge1 = [(1:size(node1,1))',[(2:size(node1,1))'; 1]];
      edge2 = [1,2; 2,3; 3,4; 4,1];

      edge = [edge1; edge2+size(node1,1)];
      node = [node1; node2];

      face{1} = 1:size(edge1,1);
      face{2} = 1:size(edge,1);
      
      hdata.hmax = 1.0;

      [p,t]=meshfaces(node,edge,face,hdata);
      
    case 3
        % Geometry
        dtheta = pi/6; theta = (-pi:dtheta:(pi-dtheta))';
        R = 0.2;
        x1 = 0.25; y1 = 0.25;
        x2 = 0.75; y2 = 0.25;
        x3 = 0.75; y3 = 0.75;
        x4 = 0.25; y4 = 0.75;
        node1 = [x1 + R*cos(theta), y1 + R*sin(theta)];
        node2 = [x2 + R*cos(theta), y2 + R*sin(theta)];
        node3 = [x3 + R*cos(theta), y3 + R*sin(theta)];
        node4 = [x4 + R*cos(theta), y4 + R*sin(theta)];
        nodeb = [0,0;1,0;1,1;0,1];
        
        n1 = size(node1,1);
        n2 = size(node2,1);
        n3 = size(node3,1);
        n4 = size(node4,1);
        
        edge1 = [(1:size(node1,1))',[(2:size(node1,1))'; 1]];
        edge2 = [n1+(1:size(node2,1))',[n1+(2:size(node2,1))'; 1+n1]];
        edge3 = [n1+n2+(1:size(node3,1))',[n1+n2+(2:size(node3,1))'; 1+n1+n2]];
        edge4 = [n1+n2+n3+(1:size(node4,1))',[n1+n2+n3+(2:size(node4,1))'; 1+n1+n2+n3]];
        edgeb = [n1+n2+n3+n4+(1:size(nodeb,1))',[n1+n2+n3+n4+(2:size(nodeb,1))'; 1+n1+n2+n3+n4]];
        
%         edge = [edge1;edge3;edgeb];
%         node = [node1;node3;nodeb];
        edge = [edge1;edge2;edge3;edge4;edgeb];
        node = [node1;node2;node3;node4;nodeb];
        
%         face{1} = 1:size(edge1,1);
%         face{2} = 1:size(edge3,1)+size(edge1,1);
%         face{3} = 1:size(edge,1);
        face{1} = 1:size(edge1,1);
        face{2} = 1:size(edge2,1)+size(edge1,1);
        face{3} = 1:size(edge3,1)+size(edge1,1)+size(edge2,1);
        face{4} = 1:size(edge4,1)+size(edge1,1)+size(edge2,1)+size(edge3,1);
        face{5} = 1:size(edge,1);
        
        options.mlim = 0.001;
        
        [p,t,fnum]=meshfaces(node,edge,face);

   otherwise
      error('Invalid demo. N must be between 1-3');
end

end      % facetest()
