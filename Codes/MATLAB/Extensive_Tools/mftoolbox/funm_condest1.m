function [c,est] = funm_condest1(A,fun,fun_frechet,flag1,varargin)
%FUNM_CONDEST1  Estimate of 1-norm condition number of matrix function.
%    C = FUNM_CONDEST1(A,FUN,FUN_FRECHET,FLAG) produces an estimate of
%    the 1-norm relative condition number of function FUN at the matrix A.
%    FUN and FUN_FRECHET are function handles:
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
%    MATLAB'S NORMEST1 (block 1-norm power method) is used, with a random
%    starting matrix, so the approximation can be different each time.
%    C = FUNM_CONDEST1(A,FUN,FUN_FRECHET,FLAG,P1,P2,...) passes extra inputs
%    P1,P2,... to FUN and FUN_FRECHET.
%    [C,EST] = FUNM_CONDEST1(A,...) also returns an estimate EST of the
%    1-norm of the Frechet derivative.
%    Note: this function makes an assumption on the adjoint of the
%    Frechet derivative that, for f having a power series expansion,
%    is equivalent to the series having real coefficients.

if nargin < 3 || isempty(fun_frechet), fte_diff = 1; else fte_diff = 0; end
if nargin < 4 || isempty(flag1), flag1 = 0; end

n = length(A);
funA = feval(fun,A,varargin{:});
if fte_diff, d = sqrt( eps*norm(funA,1) ); end

factor = norm(A,1)/norm(funA,1);

[est,v,w,iter] = normest1(@afun);
c = est*factor;

       %%%%%%%%%%%%%%%%%%%%%%%%% Nested function.
       function Z = afun(flag,X)
       %AFUN  Function to evaluate matrix products needed by NORMEST1.

       if isequal(flag,'dim')
          Z = n^2;
       elseif isequal(flag,'real')
          Z = isreal(A);
       else

          [p,q] = size(X);
          if p ~= n^2, error('Dimension mismatch'), end
          E = zeros(n);
          Z = zeros(n^2,q);
          for j=1:q

              E(:) = X(:,j);

              if isequal(flag,'notransp')

                 if fte_diff
                    Y = (feval(fun,A+d*E/norm(E,1),varargin{:}) - funA)/d;
                 else
                    if flag1
                       Y = feval(fun_frechet,'notransp',E,varargin{:});
                    else
                       Y = feval(fun_frechet,A,E,varargin{:});
                    end
                 end

              elseif isequal(flag,'transp')

                 if fte_diff
                    Y = (feval(fun,A'+d*E/norm(E,1),varargin{:}) - funA')/d;
                 else
                    if flag1
                       Y = feval(fun_frechet,'transp',E,varargin{:});
                    else
                       Y = feval(fun_frechet,A',E,varargin{:});
                    end
                 end

              end

              Z(:,j) = Y(:);
          end

       end

       end

end
