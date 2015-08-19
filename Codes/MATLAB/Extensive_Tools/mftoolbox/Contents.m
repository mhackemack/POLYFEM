% Matrix Function Toolbox (MFT).
% Version 1.0             6-Mar-2008
% Copyright (c) 2008 by N. J. Higham
%
%   arnoldi               - Arnoldi iteration
%   ascent_seq            - Ascent sequence for square (singular) matrix.
%   cosm                  - Matrix cosine by double angle algorithm.
%   cosm_pade             - Evaluate Pade approximation to the matrix cosine.
%   cosmsinm              - Matrix cosine and sine by double angle algorithm.
%   cosmsinm_pade         - Evaluate Pade approximations to matrix cosine and sine.
%   expm_cond             - Relative condition number of matrix exponential.
%   expm_frechet_pade     - Frechet derivative of matrix exponential via Pade approx.
%   expm_frechet_quad     - Frechet derivative of matrix exponential via quadrature.
%   fab_arnoldi           - f(A)*b approximated by Arnoldi method.
%   funm_condest1         - Estimate of 1-norm condition number of matrix function.
%   funm_condest_fro      - Estimate of Frobenius norm condition number of matrix function.
%   funm_ev               - Evaluate general matrix function via eigensystem.
%   funm_simple           - Simplified Schur-Parlett method for function of a matrix.
%   logm_cond             - Relative condition number of matrix logarithm.
%   logm_frechet_pade     - Frechet derivative of matrix logarithm via Pade approx.
%   logm_iss              - Matrix logarithm by inverse scaling and squaring method.
%   logm_pade_pf          - Evaluate Pade approximant to matrix log by partial fractions.
%   mft_test              - Test the Matrix Function Toolbox.
%   mft_tolerance         - Convergence tolerance for matrix iterations.
%   polar_newton          - Polar decomposition by scaled Newton iteration.
%   polar_svd             - Canonical polar decomposition via singular value decomposition.
%   polyvalm_ps           - Evaluate polynomial at matrix argument by Paterson-Stockmeyer alg.
%   power_binary          - Power of matrix by binary powering (repeated squaring).
%   quasitriang_struct    - Block structure of upper quasitriangular matrix.
%   readme                - Welcome to the Matrix Function Toolbox.
%   riccati_xaxb          - Solve Riccati equation XAX = B in positive definite matrices.
%   rootpm_newton         - Coupled Newton iteration for matrix pth root.
%   rootpm_real           - Pth root of real matrix via real Schur form.
%   rootpm_schur_newton   - Matrix pth root by Schur-Newton method.
%   rootpm_sign           - Matrix Pth root via matrix sign function.
%   signm                 - Matrix sign decomposition.
%   signm_newton          - Matrix sign function by Newton iteration.
%   sqrtm_db              - Matrix square root by Denman-Beavers iteration.
%   sqrtm_dbp             - Matrix square root by product form of Denman-Beavers iteration.
%   sqrtm_newton          - Matrix square root by Newton iteration (unstable).
%   sqrtm_newton_full     - Matrix square root by full Newton method.
%   sqrtm_pd              - Square root of positive definite matrix via polar decomposition.
%   sqrtm_pulay           - Matrix square root by Pulay iteration.
%   sqrtm_real            - Square root of real matrix by real Schur method.
%   sqrtm_triang_min_norm - Estimated min norm square root of triangular matrix.
%   sylvsol               - Solve Sylvester equation.
