function J = get_simplex_jacobian(dim, verts)
if dim == 1
    J = verts(2) - verts(1);
else
    J = zeros(dim,dim);
    vv = verts';
    for i=1:dim
        J(:,i) = vv(:,i+1) - vv(:,1);
    end
end