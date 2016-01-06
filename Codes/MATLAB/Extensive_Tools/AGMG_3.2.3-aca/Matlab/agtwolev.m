function [P ind] = agtwolev(A,l,sym,npass,kappa,checkdd,targetcf,fracnz,trspos,verbose)
%
% Compute aggregation and associated prolongation matrix according to the
% algorithms in [2] (symmetric matrices) or [3] (general matrices)
%
% [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
%    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012.
%
% [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
%    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012.
%
% USAGE:
% [P Ind] = agtwolev(A)
% [P Ind] = agtwolev(A,l,sym,npass,kappa,checkdd,targetcoarsefac
%                   ,fracnegrcsum,trspos,verbose)
%
% INPUT:
%
%  A: sparse matrix to deal with.
%
%  l:  l==1: (default) top level aggregation
%            (priority rule according to CMK permutation)
%      l==2: further level aggregation (priority according to input ordering).
%
%  sym: sym==0: general matrix (default); sym==1: symmetric matrix
%
%  npass   is the maximal number of pairwise aggregation passes (default:2).
%
%  kappa is the threshold used to accept or not a tentative aggregate
%        (default: 10 (sym==0) or 8 (sym==1)).
%
%  checkdd is the threshold to keep outside aggregation nodes where
%         the matrix is strongly diagonally dominant (based on mean of row
%         and column);
%         In fact, uses the maximum of |checkdd| and of the default value
%            according to kappa as indicated in [1,2]
%            (hence |checkdd| < 1 ensures that one uses this default value)
%         checkdd <0 : consider |checkdd|, but base the test on minus the
%               sum of offdiagonal elements, without taking the absolute value
%         (default: 0.5).
%
%  targetcoarsefac is the target coarsening factor (parameter tau in the main
%         coarsening algorithm in [1,2]): further pairwise aggregation passes
%         are omitted once the number of nonzero entries has been reduced by a
%         factor of at least targetcoarsefac (default: 2^npass).
%
%  fracnegrcsum: if, at some level, more than fracnegrcsum*nl nodes,
%         where nl is the total number of nodes at that level, have
%         negative mean row and column sum, then the aggregation algorithm
%         of [2,3] is modified, exchanging all diagonal entries for the mean
%         row and column sum (that is, the algorithm is applied to a
%         modified matrix with mean row and column sum enforced to be zero);
%         (default: 0.25).
%
%  trspos is a threshold: if a row has a positive offdiagonal entry larger
%         than trspos times the diagonal entry, the corresponding node is
%         transferred unaggregated to the coarse grid (default: 0.45).
%
% OUTPUT:
%  ind: vector of length N, where N is the numbers of rows & columns in A;
%       ind(i) is the index of the aggregates to which i
%       belongs; ind(i)=0 iff i has been kept outside aggregation (see checkdd).
%    P: associated prolongation matrix; sparse N x Nc matrix, where Nc is the
%       number of aggregates;
%       P(i,j)=1 iff ind(i)=j and P(i,j)=0 otherwise (i=1,N ; j=1,Nc).
%
% COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB)
% This function is part of AGMG software package, release 3.2.3
% ALL USAGE OF AGMG IS SUBJECT TO LICENSE
% Enter agtwolev() for detailed conditions of use.
% See the Web pages <http://homepages.ulb.ac.be/~ynotay/agmg> for
%    release information, a detailed userguide and possible upgrades.
%
if (nargin <1)
fprintf(' COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB) \n')
fprintf(' The function agtwolev is part of AGMG software package \n')
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
fprintf(' agtwolev requires at least one input argument \n')
fprintf(' enter "help agtwolev" for more information \n')
   return
end
if  (isempty(A) || any(class(A)~='double') || ~issparse(A) || ~isreal(A))
     error('MATLAB:agmg:NonSparseMatrix', ...
    'Input Matrix A must be a sparse nonempty array of double real');
end
if nargin < 2,    l   =  1; end
if nargin < 3,    sym   =0; end;
if sym       ,    sym   =1; else, sym=0; end
if nargin < 4,    npass =2; end
if nargin < 5,    if (sym), kappa=8; else, kappa=10; end; end
if nargin < 6,    checkdd = 0.5;  end
if nargin < 7,    targetcf= 2^npass;    end
if nargin < 8,    fracnz  = 0.25; end
if nargin < 9,    trspos  = 0.45; end
if nargin < 10,    verbose = 6;    end
if verbose   ,    verbose = 6;    end
ind = dmtlagtwolev(A,l,sym,npass,kappa,checkdd,targetcf,fracnz,trspos,verbose);
%%%%%%%%%%%%%%%%%%%%%%%
ind2 = find(ind);
P = sparse(ind2, ind(ind2), sign(ind(ind2)), size(A, 1), max(ind));
%%%%%%%%%%%%%%%%%%%%%%%
return;
