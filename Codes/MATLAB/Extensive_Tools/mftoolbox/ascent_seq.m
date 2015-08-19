function [d,a] = ascent_seq(A)
%ASCENT_SEQ   Ascent sequence for square (singular) matrix.
%   [d,a] = ASCENT_SEQ(A) returns symbolically computed
%   a(i) = dim(null(A^(i-1))) and the ascent sequence d = DIFF(a).
%   A has a square root if in the ascent sequence no two terms are
%   the same odd integer.
%   This function is intended for singular matrices of small
%   dimension with exactly known entries.
%   It requires the Symbolic Math Toolbox.

if isempty(ver('symbolic'))
   error('The Symbolic Math Toolbox is required.')
end

n = length(A);
a = zeros(n,1);
A = sym(A);
X = sym(eye(n));
for i = 2:n+1
    X = X*A;
    N = null(X);
    if isempty(N), break, end
    a(i) = rank(N);
end
d = diff(a);
