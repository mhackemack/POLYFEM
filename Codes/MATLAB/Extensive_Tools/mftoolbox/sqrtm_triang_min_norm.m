function R = sqrtm_triang_min_norm(T)
%SQRTM_TRIANG_MIN_NORM  Estimated min norm square root of triangular matrix.
%   R = SQRTM_TRIANG_MIN_NORM(T) computes a primary square root of the
%   upper triangular matrix T and attempts to minimize its 1-norm.

if ~isequal(T,triu(T)), error('T must be upper triangular'), end

n = length(T);
rp = zeros(n,1);
rm = zeros(n,1);

R = zeros(n);
for j=1:n
    rp(j) = sqrt(T(j,j));
    rm(j) = -sqrt(T(j,j));
    for i=j-1:-1:1
        rp(i) = (T(i,j) - R(i,i+1:j-1)*rp(i+1:j-1))/(R(i,i) + rp(j));
        rm(i) = (T(i,j) - R(i,i+1:j-1)*rm(i+1:j-1))/(R(i,i) + rm(j));
    end
    if norm(rp(1:j),1) <= norm(rm(1:j),1)
       R(1:j,j) = rp(1:j);
    else
       R(1:j,j) = rm(1:j);
    end
end
