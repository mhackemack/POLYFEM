function mft_test(n)
%MFT_TEST  Test the Matrix Function Toolbox.
%   MFT_TEST(N) tests most of the functions in the
%   Matrix Function Toolbox on (mainly) random matrices of order N.
%   The default is N = 8.
%   For a run time of a few seconds choose a small
%   value of N (such as the default).
%   Each invocation uses different random matrices.
%   This is not a thorough suite of tests.

%   For N much larger than 10 it may be necssary to adjust the fudge factors
%   in this function and in MFT_TOLERANCE in order to achieve successful
%   completion.

if nargin < 1, n = 8; end

fudge_factor = 200;
A = randn(n);
E = randn(n);
B = A^2;  % No eigenvalues on nonpositive real axis.
m = 2*n;
c = randn(m,1);
Ap = A'*A; Bp = A'*A; % Positive definite.

C = randn(n);
X = sylvsol(A,E,C);
assert_small( (A*X+X*E-C)/((norm(A,1)+norm(E,1))*norm(X,1) + norm(C,1)))
assert_eq(real(X),X)

T = triu(randn(n)); U = triu(randn(n));
X = sylvsol(T,U,C);
assert_eq(real(X),X)
assert_small( (T*X+X*U-C)/((norm(T,1)+norm(U,1))*norm(X,1) + norm(C,1)))

if license('test','Symbolic_Toolbox') && ~isempty(ver('symbolic'))
   % If Symbolic Math Toolbox licensed and installed...
    for p = 2:6
      J = gallery('jordbloc',p,0);
      [d,a] = ascent_seq(J);
      assert_eq(d,ones(p,1))
      assert_eq(a,(0:p)')
    end
end

A1 = 2.5*A/sqrt(norm(A^2,inf));
for k = 1:3  % Try different norms.
  C1 = cosm(A1);
  [C,S] = cosmsinm(A1);
  C0 = funm(A1,@cos);
  S0 = funm(A1,@sin);
  assert_sim(C,C1)
  assert_sim(C,C0)
  assert_sim(S,S0)
  assert_small( (C^2+S^2-eye(n)) / (norm(C,1)^1+norm(S,1)^2) )
  A1 = 2*A1;
end

L = expm_frechet_pade(A,E);
R = expm_frechet_quad(A,E);
assert_sim(L,R,1e-2)

% Check identity $L_{\exp}\bigl(\log(A), L_{\log}(A,E)\bigr) = E$.
L_log = logm_frechet_pade(B,E);
L_exp = expm_frechet_pade(logm(B),L_log);
assert_sim(L_exp,E,sqrt(eps)*norm(E,1))

m = n-1;
E = eye(m);
E = E(:,m);
[Q,H] = arnoldi(A,randn(n,1),m);
res = A*Q(:,1:m) - Q(:,1:m)*H(1:m,1:m) - H(m+1,m)*Q(:,m+1)*E';
assert_small(res/(norm(H,1)*norm(Q,1)))

b = randn(n,1); y = fab_arnoldi(A,b,@exp,n);
assert_sim(y,expm(A)*b)

X = logm(B);
X1 = logm_iss(B);
tol = funm_condest1(B,@logm,@logm_frechet_pade)*eps*n;
assert_sim(X,X1,tol)

c1 = expm_cond(A);
[c1a,K] = expm_cond(A);
assert_sim(c1,c1a)
st = randn('state');
c2 = funm_condest_fro(A,@expm,@expm_frechet_pade);
randn('state',st);
c3 = funm_condest_fro(A,@expm,@fun_frechet_exp,[],1);
assert_sim(c2,c3,0.5)
fudge_factor1 = 2;
if c2 < c1/10 || c2 > c1*fudge_factor1, [c2 c1 c2/c1], error('Failure'), end
if c3 < c1/10 || c3 > c1*fudge_factor1, [c3 c1 c3/c1], error('Failure'), end

c1 = funm_condest1(A,@expm);
st = rand('twister');
c2 = funm_condest1(A,@expm,@expm_frechet_pade);
assert_sim(c1,c2,1)
rand('twister',st);
c3 = funm_condest1(A,@expm,@fun_frechet_exp,1);
assert_sim(c2,c3,0.5)

c1 = funm_condest_fro(B,@logm);
c2 = funm_condest_fro(B,@logm,@logm_frechet_pade);
assert_sim(c1,c2,1)

[U1,H1] = polar_newton(A);
assert_small((A-U1*H1)/norm(A,1))
[U2,H2] = polar_svd(A);
assert_small((A-U2*H2)/norm(A,1))
assert_sim(U1,U2)
assert_small((H1-H2)/norm(H1,1))
assert_small(U1'*U1-eye(n))
assert_small(U2'*U2-eye(n))
assert_eq(H1,H1')
assert_eq(H2,H2')

A2 = randn(n+2,n);
for i = 1:4
   [U,H] = polar_svd(A2);
   assert_small((A2-U*H)/norm(A2,1))
   assert_small(U1'*U1-eye(n))
   assert_eq(H,H')
   A2 = A2';
   if i == 2, A2(:,round(n/2)) = 0; end
end

P1 = polyvalm_ps(c,A);
P2 = polyvalm(c,A);
assert_small((P1-P2)/norm(P1,1))

assert_sim(power_binary(A,m),A^m)

X = riccati_xaxb(Ap,Bp);
assert_small( (X*Ap*X-Bp) / (norm(X,1)^2*norm(Ap,1)+norm(Bp,1)) )
assert_eq(X,X')

for p = [2 5 10 16]
  X = rootpm_real(B,p); assert_small( (X^p-B)/(norm(X,1)^p + norm(B,1)) );
  X = rootpm_sign(B,p); assert_small( (X^p-B)/(norm(X,1)^p + norm(B,1)) );
  [X,Y] = rootpm_schur_newton(B,p);
  assert_small( (X^p-B)/(norm(X,1)^p + norm(B,1)) )
  assert_small( (X*Y - eye(n))/(norm(X,1)*norm(Y,1)) )
end

[S,N] = signm(A);
assert_small( (S^2 - eye(n))/(norm(S,1)^2+1) )
assert_small( (A-S*N)/(norm(A,1)+norm(S,1)*norm(N,1)) )

[X0,alpha,condest] = sqrtm(B);
tol = n*norm(X0,1)*condest*eps;

[P,Q] = sqrtm_db(B);
assert_small( (P^2 - B)/(norm(P,1)^2+norm(B,1)) )
assert_small( (Q^2 - inv(B))/(norm(Q,1)^2+norm(inv(B),1)) )
assert_small(X0-P, tol)
[P,Q] = sqrtm_dbp(B);
assert_small( (P^2 - B)/(norm(P,1)^2+norm(B,1)) )
assert_small( Q - eye(n) )
assert_small(X0-P, tol)

% Following usually succeeds but can fail:
% only local cgce conditions are known for full Newton.
% X = sqrtm_newton_full(B);
% assert_small( (X^2 - B)/(norm(X,1)^2+norm(B,1)) )

X = sqrtm_pd(Ap);
assert_small( (X^2 - Ap)/(norm(X,1)^2+norm(Ap,1)) )
assert_eq(X,X')

C = full(gallery('tridiag',n,1,4,1));
D = diag(diag(C));
X = sqrtm_pulay(C,D);
assert_small( (X^2 - C)/(norm(X,1)^2+norm(C,1)) )

X = sqrtm_real(B);
assert_small( (X^2 - B)/(norm(X,1)^2+norm(B,1)) )
assert_small(X0-P, tol)

T = schur(A,'complex');
R = sqrtm_triang_min_norm(T);
assert_small( (R^2 - T)/(norm(R,1)^2+norm(T,1)) )

fprintf(['MFT_TEST: All tests of the Matrix Function Toolbox passed' ...
        ' (n = %g).\n'], n)

      % Nested functions

      function L = fun_frechet_exp(flag,E)
      % Frechet derivative of exponential.

      if strcmp(flag,'transp'), E = E'; end
      L = expm_frechet_pade(A,E);
      if strcmp(flag,'transp'), L = L'; end

      end

      % ---------------------------------------------------------
      % Assertion functions.

      function assert_sim(a,b,tol)
      if nargin < 3, tol = fudge_factor*eps(superiorfloat(a,b))*length(a); end
      if norm(a-b,1)/max( norm(a,1), norm(b,1) ) > tol
         fprintf('%9.2e, %9.2e\n', ...
                  norm(a-b,1)/max( norm(a,1), norm(b,1)), tol )
         error('Failure')
      end
      end

      function assert_small(a,tol)
      if nargin < 2, tol = fudge_factor*eps(class(a))*length(a); end
      if norm(a,1) > tol
         fprintf('%9.2e, %9.2e\n',norm(a,1), tol),
         error('Failure')
      end
      end

      function assert_eq(a,b)
      if norm(a-b,1)
         error('Failure')
      end
      end

end
