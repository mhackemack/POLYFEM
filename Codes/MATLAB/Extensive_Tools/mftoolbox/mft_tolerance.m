function tol = mft_tolerance(A)
%MFT_TOLERANCE   Convergence tolerance for matrix iterations.
%   TOL = MFT_TOLERANCE(A) returns a convergence tolerance to use in
%   the matrix iterations in the Matrix Function Toolbox applied to the
%   matrix A.  All functions in the toolbox call this function to set
%   the convergence tolerance.

tol = sqrt(length(A))*eps/2;
