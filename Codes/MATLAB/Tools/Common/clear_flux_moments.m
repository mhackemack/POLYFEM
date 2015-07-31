function x = clear_flux_moments(x, n)
[mx,nx] = size(x);
for i=1:mx
    for j=1:nx
        x{i,j} = zeros(n,1);
    end
end