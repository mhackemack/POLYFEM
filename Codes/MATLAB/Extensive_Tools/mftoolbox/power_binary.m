function X = power_binary(A,m)
%POWER_BINARY   Power of matrix by binary powering (repeated squaring).
%   X = POWER_BINARY(A,m) computes A^m for a square matrix A and a
%   positive integer m, by binary powering.

s = double(dec2bin(m)) - 48;  % Binary representation of s in double array.
k = length(s);

P = A;
i = k;
while s(i) == 0
      P = P^2;
      i = i-1;
end
X = P;
for j = i-1:-1:1
    P = P^2;
    if s(j) == 1
       X = X*P;
    end

end
