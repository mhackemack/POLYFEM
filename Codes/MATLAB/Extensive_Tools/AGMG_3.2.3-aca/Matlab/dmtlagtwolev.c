/*
! COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB)
!
! This file is part of AGMG software package
! Release 3.2.3-aca built on "Dec  8 2015" by Yvan Notay
!
! ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE".
! IF YOU OBTAINED A COPY OF THIS SOFTWARE WITHOUT THIS FILE,
! PLEASE CONTACT ynotay@ulb.ac.be
!
! In particular, if you have a free academic license:
!
! (1) You must be a member of an educational, academic or research institution.
!     The license agreement automatically terminates once you no longer fulfill
!     this requirement.
!
! (2) You are obliged to cite AGMG in any publication or report as:
!     "Yvan Notay, AGMG software and documentation;
!      see http://homepages.ulb.ac.be/~ynotay/AGMG".
!
! (3) You may not make available to others the software in any form, either
!     as source or as a precompiled object.
!
! (4) You may not use AGMG for the benefit of any third party or for any
!     commercial purposes. Note that this excludes the use within the
!     framework of a contract with an industrial partner.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! See the Web pages <http://homepages.ulb.ac.be/~ynotay/AGMG> for
! release information and possible upgrade.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DICLAIMER:
!    AGMG is provided on an "AS IS" basis, without any explicit or implied
!    WARRANTY; see the file "LICENSE" for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   If you use AGMG for research, please observe that your work benefits
!   our past research efforts that allowed the development of AGMG.
!   Hence, even if you do not see it directly, the results obtained thanks
!   to the use of AGMG depend on the results in publications [1-3] below,
!   where the main algorithms used in AGMG are presented and justified.
!   It is then a normal duty to cite these publications (besides citing
!   AGMG itself) in any scientific work depending on the usage of AGMG,
!   as you would do with any former research result you are using.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
#include "mex.h"
#include "matrix.h"

/* computational subroutines */
extern void dag2l_twolev_(int *, int *, double *, int *, int *, int *,
			   int *, int *, double *, double *, double *,
			   double *, double *, int *);


/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    mwIndex *irs,*jcs;
    double  *as, *kap, *checkd, *targetcf, *fracnz, *trspos;
    double  *npas1, *sym1, *l1, *x , *iprint1;
    mwSize  nn, nsiz, ncol;
    int     *ja, *ia, *ind, n, nz, i, l, sym, npas, iprint;
  
    /* get pointers to input*/
    nn      = mxGetN(prhs[0]);
    as      = mxGetPr(prhs[0]);
    irs     = mxGetIr(prhs[0]);
    jcs     = mxGetJc(prhs[0]);
    l1      = mxGetPr(prhs[1]);
    sym1    = mxGetPr(prhs[2]);
    npas1   = mxGetPr(prhs[3]);
    kap     = mxGetPr(prhs[4]);
    checkd  = mxGetPr(prhs[5]);
    targetcf= mxGetPr(prhs[6]);
    fracnz  = mxGetPr(prhs[7]);
    trspos  = mxGetPr(prhs[8]);
    iprint1 = mxGetPr(prhs[9]);
    n       = nn;
    
    /* transform to integer */
    l=(int)*l1;
    sym=(int)*sym1;
    npas=(int)*npas1;
    iprint=(int)*iprint1;
  
    /* create workspace, get pointers to it*/
    nz = jcs[n];
    nsiz = nz*sizeof(int);
    ja  = mxMalloc(nsiz);
    nsiz = (n+1)*sizeof(int);
    ia  = mxMalloc(nsiz);
        for (i=0 ; i<=n ; i++) {
            ia[i] = jcs[i]+1;}
        for (i=0 ; i<jcs[n] ; i++) {
            ja[i] = irs[i]+1;}
    nsiz = n*sizeof(int);
    ind    = mxMalloc(nsiz);
    /* call the aggregation procedure */
    dag2l_twolev_(&l,&n,as,ja,ia,ind,&sym,&npas,kap,checkd,targetcf,fracnz,trspos,&iprint);

    /* create a new array and set the output pointer to it */
    ncol=1;
    plhs[0] = mxCreateDoubleMatrix(nn, ncol, mxREAL);
    x = mxGetPr(plhs[0]);
        for (i=0 ; i<n ; i++) {
            x[i] = ind[i];
        }
}
