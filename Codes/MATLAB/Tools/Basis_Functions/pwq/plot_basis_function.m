function plot_basis_function(fn,X,Y,sides)

[n,m]=size(X);
Z=zeros(n,m);
for i=1:n
    for j=1:m
        Z(i,j)=fn(X(i,j),Y(i,j),sides);
    end
end

surf(Y,X,Z');
xlabel('X');
ylabel('Y');
view(0,90)