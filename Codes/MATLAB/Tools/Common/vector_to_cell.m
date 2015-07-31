function out = vector_to_cell(x,DoF)
ndof = DoF.TotalDoFs;
out = cell(length(x)/ndof,1);
for i=1:length(out)
    out{i} = x((i-1)*ndof+1:i*ndof);
end