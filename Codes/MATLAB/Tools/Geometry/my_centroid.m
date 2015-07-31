function [C] = my_centroid( verts )
dim = size(verts,2);
if dim == 2
    x=verts(:,1);
    y=verts(:,2);
    aux = x.*y([2:end 1]) - y.*x([2:end 1]);
    signed_area=0.5* sum(aux);
    C(1) = sum( ( x + x([2:end 1]) ).*aux ) / (6*signed_area);
    C(2) = sum( ( y + y([2:end 1]) ).*aux ) / (6*signed_area);
elseif dim == 3
    
end