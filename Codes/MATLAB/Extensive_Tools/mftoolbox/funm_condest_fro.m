function c = funm_condest_fro(A,fun,fun_frechet,its,flag,varargin)
%FUNM_CONDEST_FRO  Estimate of Frobenius norm condition number of matrix function.
%    C = FUNM_CONDEST_FRO(A,FUN,FUN_FRECHET,ITS,FLAG) produces an estimate of
%    the Frobenius norm relative condition number of function FUN at
%    the matrix A.    FUN and FUN_FRECHET are function handles:
%      - FUN(A) evaluates the function at the matrix A.
%      - If FLAG == 0 (default)
%           FUN_FRECHET(B,E) evaluates the Frechet derivative at B
%              in the direction E;
%        if FLAG == 1
%           - FUN_FRECHET('notransp',E) evaluates the
%                    Frechet derivative at A in the direction E.
%           - FUN_FRECHET('transp',E) evaluates the
%                    Frechet derivative at A' in the direction E.
%    If FUN_FRECHET is empty then the Frechet derivative is approximated
%    by finite differences.  More reliable results are obtained when
%    FUN_FRECHET is supplied.
%    The power method is used, with a random starting matrix,
%    so the approximation can be different each time.
%    ITS iterations (default 6) are used.
%    C = FUNM_COND_EST_FRO(A,FUN,FUN_FRECHET,ITS,FLAG,P1,P2,...)
%    passes extra inputs P1,P2,... to FUN and FUN_FRECHET.
%    Note: this function makes an assumption on the adjoint of the
%    Frechet derivative that, for f having a power series expansion,
%    is equivalent to the series having real coefficients.

if nargin < 5 || isempty(flag), flag = 0; end
if nargin < 4 || isempty(its), its = 6; end
if nargin < 3 || isempty(fun_frechet), fte_diff = 1; else fte_diff = 0; end

funA = feval(fun,A,varargin{:});
d = sqrt( eps*norm(funA,'fro') );
Z = randn(size(A));
Znorm = 1;

factor = norm(A,'fro')/norm(funA,'fro');

for i=1:its

   Z = Z/norm(Z,'fro');
   if fte_diff
      W = (feval(fun,A+d*Z,varargin{:}) - funA)/d;
   else
      if flag
         W = feval(fun_frechet,'notransp',Z,varargin{:});
      else
         W = feval(fun_frechet,A,Z,varargin{:});
      end
   end

   W = W/norm(W,'fro');
   if fte_diff
      Z = (feval(fun,A'+d*W,varargin{:}) - funA')/d;
   else
      if flag
         Z = feval(fun_frechet,'transp',W,varargin{:});
      else
         Z = feval(fun_frechet,A',W,varargin{:});
      end
   end

   Znorm = norm(Z,'fro');

end

c = Znorm*factor;
