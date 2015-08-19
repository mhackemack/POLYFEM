function F = funm_ev(A,fun)
%FUNM_EV   Evaluate general matrix function via eigensystem.
%   F = FUNM_EV(A,FUN) evaluates the function FUN at the
%   square matrix A using the eigensystem of A.
%   This function is intended for diagonalizable matrices only
%   and can be numerically unstable.

[V,D] = eig(A);
F = V * diag(feval(fun,diag(D))) / V;
