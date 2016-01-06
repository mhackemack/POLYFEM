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
extern void zagmg_(int *,double *,int *,int *,double *,double *,int *,int *,
                   int *,int *,double *);

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    mwIndex *irs,*jcs;
    double  *as, *ai, *xr, *xi, *x, *fr, *fi, *fc;
    double  *tolp, *iprintp, *iterp, *nrestp, *x0r, *x0i, *ijobp;
    double  tol;
    mwSize  nn, ncol, nsiz;
    int     n, i, nrest, iter, iprint, ijob, ijb, ijbe;
    static int    *ja, *ia, nz, np=0;
    static double *a; 
  
    /* get pointers to input*/
    fr      = mxGetPr(prhs[1]);
    fi      = mxGetPi(prhs[1]);
    iprintp = mxGetPr(prhs[2]);
    nrestp  = mxGetPr(prhs[3]);
    iterp   = mxGetPr(prhs[4]);
    tolp    = mxGetPr(prhs[5]);
    ijobp   = mxGetPr(prhs[6]);
    iprint=(int)*iprintp;
    nrest =(int)*nrestp;
    iter  =(int)*iterp;
    tol   =*tolp;

    ijb  =(int)*ijobp;
        if (ijb >=100) {ijob=ijb-100;ijbe=ijob;}
    else {ijob=ijb;
        if (ijb >=0) {ijbe=ijb+100;}
        else {ijbe=ijb;}
    }
  
    if (ijob <= 1 && np>0) {
        np=0;
    }
    if (ijob >= 0) {
  
    if (ijob < 3) {
    nn      = mxGetN(prhs[0]);
    as      = mxGetPr(prhs[0]);
    ai      = mxGetPi(prhs[0]);
    irs     = mxGetIr(prhs[0]);
    jcs     = mxGetJc(prhs[0]);
    /* create workspace, get pointers to it*/
    n     = nn;
    nz  = jcs[n];
    nsiz= 2*nz*sizeof(double);
    a   = mxMalloc(nsiz); 
    nsiz= nz*sizeof(int);
    ja  = mxMalloc(nsiz); 
    nsiz= (n+1)*sizeof(int);
    ia  = mxMalloc(nsiz); 
    for (i=0 ; i<=n ; i++) {
          ia[i] = jcs[i]+1;}
    for (i=0 ; i<jcs[n] ; i++) {
            ja[i] = irs[i]+1;
 	    a[2*i] = as[i];
	    a[2*i+1] = ai[i];}
    }
    else {n=np;nn=np;}

    if (ijob != 1) {
    /* create workspace, get pointers to it*/
    nsiz= 2*n*sizeof(double);
    fc  = mxMalloc(nsiz);
    x   = mxMalloc(nsiz);
    /* copy r.h.s. to avoid overwriting */
    for(i=0 ; i<nn ; i++){
       fc[2*i]=fr[i];
       if (fi == NULL) {fc[2*i+1]=0;} else {fc[2*i+1]=fi[i];}
    }
    /* process initial guess, if any */
       if (nrhs > 7 && ijob < 3) {
        x0r = mxGetPr(prhs[7]);
        x0i = mxGetPi(prhs[7]);
	    for(i=0 ; i<nn ; i++){
		x[2*i]=x0r[i];
                if (x0i == NULL) {x[2*i+1]=0;} else {x[2*i+1]=x0i[i];}
	    }
	    ijbe=ijbe+10;
       }
    }
    }

         zagmg_(&n,a,ja,ia,fc,x,&ijbe,&iprint,&nrest,&iter,&tol); 

    if (ijob == 0 || ijob >1) {
       if (ijob == 3) {iter=0;}
       /* create a new array and set the output pointer to it */
       ncol = 1;
       plhs[0] = mxCreateDoubleMatrix(nn, ncol, mxCOMPLEX);
       xr = mxGetPr(plhs[0]);
       xi = mxGetPi(plhs[0]);
       for(i=0 ; i<nn ; i++){
		xr[i]=x[2*i];
                xi[i]=x[2*i+1];
	    }
       plhs[1] = mxCreateDoubleMatrix(ncol, ncol, mxREAL);
       iterp   = mxGetPr(plhs[1]);
       *iterp  = iter;
       if (iter < 0) {iter=-iter;}
       nn = iter+1;
       plhs[2] = mxCreateDoubleMatrix(nn, ncol, mxREAL);
       fr       = mxGetPr(plhs[2]);
       for (i=0 ; i <= iter; i++) {fr[i]=fc[2*i];}
    }
    /* keep n in local memory if ijob == 1 */
    else if (ijob == 1) {
        np=n;
        }
}




