function [x,flag,relres,iter,resvec]=agmg(A,b,icg,tol,maxit,verbose,x0,ijob)
%AGMG: Aggregation-based Multigrid iterative Method.
%   X = AGMG(A,B) attempts to solve the system of linear equations A*X=B for
%   X. The N-by-N coefficient matrix A must be square and the right
%   hand side column vector B must have length N.
%   The matrix A may be real or complex; for complex matrices, however,
%   AGMG is tentative only. If A is real, B has to be real too.
%
%   X = AGMG(A,B,RESTART)
%       If RESTART==1, use CG and performs some simplifications based on
%             the assumption that thecoefficient matrix A is symmetric;
%             should be used only when A is symmetric and positive definite.
%       If RESTART>=2, use GCR restarted each RESTART iterations.
%       If RESTART==0 or [] then AGMG use the default, which is
%                     GCR restarted each 10 iterations.
%
%   X = AGMG(A,B,RESTART,TOL) specifies the tolerance of the method.
%       If TOL is [] then AGMG uses the default, 1e-6.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT) specifies the maximum number of iterations.
%       If MAXIT is [] then AGMG uses the default, 100.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE)
%       If VERBOSE==1, information is displayed on the solution process
%       If VERBOSE==0 or [], AGMG works silently.
%       See the README file provided with AGMG for a description of the
%       verbose output.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE,X0) specifies an initial guess.
%       If X0 is [] then AGMG uses the default, an all zero vector.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB)
%       If IJOB==1, performs the setup only (preprocessing:
%                   prepares all parameters for subsequent solves).
%          Then, only A and VERBOSE input parameters are significant
%                and the calling may be: AGMG(A,[],[],[],[],[],[],1)
%          The returned X is empty and other output parameters are meaningless.
%       If IJOB==2, solves only, based on previous setup.
%          Then, A may differ from the matrix supplied for set up (former
%                call with IJOB==1), but it means using a preconditioner
%                computed for a matrix to solve a system with another matrix,
%                which is not recommended in general.
%       If IJOB==3, the returned X is not the solution of the linear system,
%                   but the result of the action of the multigrid
%                   preconditioner on the right hand side B.
%          Then, A, TOL, MAXIT, and X0 are not significant
%                and RESTART determines only the type of inner iterations;
%                the calling may be: AGMG([],B,[],[],[],[],[],3)
%          Further output parameters (besides X) are meaningless.
%       If IJOB==-1, erases the setup and releases internal memory.
%          Other input parameters are not significant and the calling may be
%                AGMG(A,[],[],[],[],[],[],-1)
%          The returned X is empty and other output parameters are meaningless.
%       If IJOB==0 or [], performs setup + solve + memory release
%          (default: other input & ouput parameters have their usual meaning).
%
%       IJOB == 100,110,101,102,112: same as, respectively, IJOB==0,10,1,2,12
%       but, use the TRANSPOSE of the input matrix in A, JA, IA.
%       Hence, AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB+100)
%       is equivavelent to AGMG(A',B,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB),
%       but less memory consuming and often faster.
%
%   !!! IJOB==2,3,12,102,112 require that one has previously called AGMG
%       with IJOB==1 or IJOB==101
%
%   [X,FLAG] = AGMG(A,B,...) also returns a convergence FLAG:
%    0 AGMG converged to the desired tolerance TOL within MAXIT iterations
%    1 AGMG iterated MAXIT times but did not converge.
%
%   [X,FLAG,RELRES] = AGMG(A,B,...) also returns the relative residual
%    NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%    It may happen that RELRES > TOL while FLAG is zero when TOL is
%    below the accuracy that can be attained with a backward stable solver.
%    Then RELRES is comparable to what can be obtained with a direct solver.
%    (RELRES cannot be computed when IJOB==2 since the first argument may
%     be different from the system matrix).
%
%   [X,FLAG,RELRES,ITER] = AGMG(A,B,...) also returns the iteration number
%    at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = AGMG(A,B,...) also returns a vector of the
%    residual norms at each iteration including NORM(B-A*X0).
%
%   Example: m=50; A=delsq(numgrid('S',m+2)); n=m*m; b=ones(n,1);
%            x=agmg(A,b,1);         % A is symmetric positive definite
%            x=agmg(A,b);           % Consider A as a general matrix
%            x=agmg(A,b,1,[],[],1); % Verbose output
%            % Example of tolerance below attainable accuracy:
%            [x,flag,relres,iter,resvec]=agmg(A,b,1,1e-20,[],1);
%            % AGMG report normal convergence:
%            disp('Convergence flag = '),disp(flag)
%            % The relative residual is nevertheless larger than TOL
%            % (and than reported by the verbose output):
%            disp('Relative residual = '),disp(relres)
%            % But the attained accuracy is similar to the one obtained with "\":
%            y=A\b;
%            disp('Relative residual with "\" = '),disp(norm(A*y-b)/norm(b))
%
% COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB)
% This function is part of AGMG software package, release 3.2.3
% ALL USAGE OF AGMG IS SUBJECT TO LICENSE
% Enter agmg() for detailed conditions of use.
% See the Web pages <http://homepages.ulb.ac.be/~ynotay/agmg> for
%    release information, a detailed userguide and possible upgrades.
%
persistent preprocessed n notcpl
if isempty(preprocessed)
   preprocessed = 0;
end
   % Check matrix and right hand side vector inputs have appropriate size
   if (nargin < 2)
fprintf(' COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB) \n')
fprintf(' The function agmg is part of AGMG software package \n')
fprintf(' \n')
fprintf(' ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE". \n')
fprintf(' IF YOU OBTAINED A COPY OF THIS SOFTWARE WITHOUT THIS FILE, \n')
fprintf(' PLEASE CONTACT ynotay@ulb.ac.be \n')
fprintf(' \n')
fprintf(' In particular, if you have a free academic license: \n')
fprintf(' \n')
fprintf(' (1) You must be a member of an educational, academic or research institution. \n')
fprintf('     The license agreement automatically terminates once you no longer fulfill \n')
fprintf('     this requirement. \n')
fprintf(' \n')
fprintf(' (2) You are obliged to cite AGMG in any publication or report as: \n')
fprintf('     "Yvan Notay, AGMG software and documentation; \n')
fprintf('      see http://homepages.ulb.ac.be/~ynotay/AGMG". \n')
fprintf(' \n')
fprintf(' (3) You may not make available to others the software in any form,  \n')
fprintf('     either as source or as a precompiled object. \n')
fprintf(' \n')
fprintf(' (4) You may not use AGMG for the benefit of any third party or for any \n')
fprintf('     commercial purposes. Note that this excludes the use within the \n')
fprintf('     framework of a contract with an industrial partner. \n')
fprintf(' \n')
fprintf(' See the Web pages <http://homepages.ulb.ac.be/~ynotay/AGMG> for \n')
fprintf('    release information, a detailed userguide and possible upgrades. \n')
fprintf(' \n')
fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
fprintf(' DICLAIMER: \n')
fprintf('    AGMG is provided on an "AS IS" basis, without any explicit or implied \n')
fprintf('    WARRANTY; see the see the file "LICENSE" for more details. \n')
fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
fprintf('   If you use AGMG for research, please observe that your work benefits \n')
fprintf('   our past research efforts that allowed the development of AGMG. \n')
fprintf('   Hence, even if you do not see it directly, the results obtained thanks \n')
fprintf('   to the use of AGMG depend on the results in publications [1-3] below,\n')
fprintf('   where the main algorithms used in AGMG are presented and justified.   \n')
fprintf('   It is then a normal duty to cite these publications (besides citing \n')
fprintf('   AGMG itself) in any scientific work depending on the usage of AGMG, \n')
fprintf('   as you would do with any former research result you are using. \n')
fprintf(' \n')
fprintf(' [1] Y. Notay, An aggregation-based algebraic multigrid method, \n')
fprintf('    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010 \n')
fprintf('\n')
fprintf(' [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed \n')
fprintf('    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012. \n')
fprintf(' \n')
fprintf(' [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion \n')
fprintf('    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012. \n')
fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
fprintf(' \n')
fprintf(' agmg requires at least two input arguments \n')
fprintf(' enter "help agmg" for more information \n')
   return
   end
   if (nargin < 3) || isempty(icg)
      icg = 0;
   end
   if (nargin < 4) || isempty(tol)
      tol = 1e-6;
   end
   if (nargin < 5) || isempty(maxit)
      maxit = 100;
   end
   if (nargin < 6) || isempty(verbose)
      verbose=0;
   end
   if (verbose) verbose=6; end
   if (nargin < 7)
      x0=[];
   end
   if (nargin < 8 || isempty(ijob))
      ijob=0;
   end
   if (ijob >= 100)
	   ijb=ijob-100;
   else
	   ijb=ijob;
   end
   if (ijb < -1 || ijb > 3 || round(ijb)~=ijb)
          error('MATLAB:agmg:nonvalidijob',...
'IJOB should be an integer and either IJOB or IJOB-100 should be not less than -1 and not larger than 3');
   elseif (ijb > 1 && ~preprocessed)
          error('MATLAB:agmg:nonvalidijob',...
      'Setup not done: IJOB should be equal to 0 or 100')
   elseif (ijb == -1 && ~preprocessed)
          disp('Warning: setup not done: nothing to do for IJOB == -1')
          iter=0;
          preprocessed=1;
          relres=1;
          x=[];
          resvec=[];
          return
   end
   if (ijb==0 || ijb==1 || ijb==2)
      if (isempty(A) ||  any(class(A)~='double') || ~issparse(A))
         error('MATLAB:agmg:NonSparseMatrix', ...
         'Input Matrix A must be a sparse nonempty array of double real or double complex');
      end
      if (ijb==2)
	      if (n~=size(A,1) || n~=size(A,2))
	  error('MATLAB:agmg:NotSameDim', ...
'When IJOB==2 or IJOB==102, Matrix must have same dimensions as on previous call with IJOB==1 or IJOB==101');
	      end
	      if (isreal(A) && notcpl==0)
	  error('MATLAB:agmg:NotSameType', ...
'When IJOB==2 or IJOB==102, Matrix must have same type (real or complex) as on previous call with IJOB==1 or IJOB==101');
	      end
      else
      n = size(A,1);
      if (size(A,2) ~= n)
         error('MATLAB:agmg:NonSquareMatrix', 'Matrix must be square.');
      end
      if(isreal(A)), notcpl=1; else, notcpl=0; end
      preprocessed=0;
      end
   end
   if (ijb == 0 || ijb >= 2)
      if ~isequal(size(b),[n,1])
          error('MATLAB:agmg:RSHsizeMatchCoeffMatrix', ...
          ['Right hand side must be a column vector of' ...
          ' length %d to match the coefficient matrix.'],n);
      end
      if issparse(b)
          error('MATLAB:agmg:RSHsizeMustBeFull', ...
          ['Right hand side cannot be a sparse vector.' ...
          ' Please convert it [b=full(b);] before calling agmg.']);
      end
   end
   if (notcpl)
      if (ijb == 0 || ijb >= 2)
        if(~isreal(b))
	     error('MATLAB:agmg:bcomplexwhenAreal', ...
        'For real coefficient matrix, the right hand side vector must be real');
        end
        if (isempty(x0))
           [x,iter,resvec]=dmtlagmg(A,b,verbose,icg,maxit,tol,ijob);
        else
           if(ijb == 3)
               disp('Warning: X0 ignored when IJOB==3 or IJOB==103')
           elseif(~isreal(x0))
	           disp('Warning: the coefficient matrix and the right hand side vector are real;')
               disp('         only the real part of the initial guess will be taken into account')
           end
           [x,iter,resvec]=dmtlagmg(A,b,verbose,icg,maxit,tol,ijob,x0);
         end
      else
         dmtlagmg(A,b,verbose,icg,maxit,tol,ijob);
      end
   else
      if (~preprocessed)
         disp('Warning: complex version is tentative only')
      end
      if (ijb == 0 || ijb >= 2)
        if (isempty(x0))
           [x,iter,resvec]=zmtlagmg(A,b,verbose,icg,maxit,tol,ijob);
        else
           if(ijb == 3)
               disp('Warning: X0 ignored when IJOB==3 or IJOB ==103')
           end
           [x,iter,resvec]=zmtlagmg(A,b,verbose,icg,maxit,tol,ijob,x0);
        end
      else
           zmtlagmg(A,b,verbose,icg,maxit,tol,ijob);
      end
   end
   flag=0;
   if (ijb == 1 || ijb == -1)
      iter=[];
      if (ijb == 1), preprocessed=1; else, preprocessed=0; end
      relres=[];
      x=[];
      resvec=[];
   elseif (ijb == 3)
      relres=[];
      resvec=[];
   else
      if (iter < 0), flag=1; iter=-iter; end
      if (nargout > 2), relres=norm(A*x-b)/norm(b); else, relres=[]; end
   end
end
