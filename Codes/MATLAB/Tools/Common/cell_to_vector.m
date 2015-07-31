function out = cell_to_vector(x,DoF)
n = size(x,1);
ndof = DoF.TotalDoFs;
ntot = n*ndof;
out = zeros(ntot,1);
nn = 0;
for i=1:n
    nntt = length(x{i,1});
    out(nn+1:nn+nntt) = x{i,1};
    nn = nn + nntt;
end