! COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB)
!
! This file is part of AGMG software package
! Release 3.2.3-aca built on "Dec  8 2015" by Yvan Notay
!
! ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE".
! IF YOU OBTAINED A DCOPY OF THIS SOFTWARE WITHOUT THIS FILE,
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
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! PARAMETERS DEFINITON -  may be tuned by expert users
!-----------------------------------------------------------------------
  MODULE dag2l_mem
    SAVE
!
! INTEGER
!
!  maxlev  is the maximal number of levels
!          (should be large enough - much larger than needed is armless).
!  real_len is the length of 1 REAL(kind(0.0d0)) in byte
!        (used only to display information on memory usage).
!  nsmooth  is the number of pre- and post-smoothing steps;
!  smoothtype indicates which smoother use:
!    if smoothtype==1, the smoothing scheme is
!        pre-smoothing: Forward Gauss-Seidel, then Backward GS, then Fwd GS,etc
!       post-smoothing: Backward Gauss-Seidel, then Forward GS, then Bwd GS,etc
!                       (nsmooth sweeps for both pre- and post-smoothing)
!    if smoothtype==1, ILU(0) smoother is used
!  nstep   is the maximum number of coarsening steps
!          nstep<0 means that coarsening is stopped only according to
!          the matrix size, see parameter maxcoarsesize.
!  nlvcyc  is the number of coarse levels (from bottom) on which V-cycle
!          formulation is enforced (Rmk: K-cycle always allowed at first
!          coarse level).
!  npass   is the maximal number of pairwise aggregation passes for each
!          coarsening step, according to the algorithms in [2,3].
!  maxcoarsesize: when the size of the coarse grid matrix is less than or
!                 equal to maxcoarsesize*N^(1/3),  it is factorized exactly
!                 and coarsening is stopped;
!         maxcoarsesizeslow: in case of slow coarsening,
!                 exact factorization can be performed when the size of
!                 the coarse grid matrix is less than or equal to
!                 maxcoarsesizeslow*N^(1/3).
    INTEGER, PARAMETER :: maxlev=40, real_len=8
    INTEGER, PARAMETER :: nsmooth=1, smoothtype=1, nstep=-1, nlvcyc=0
    INTEGER, PARAMETER :: maxcoarsesize=0,maxcoarsesizeslow=0
    INTEGER :: npass=2
!
! REAL
!
!  trspos is a threshold: if a row has a positive offidiagonal entry larger
!         than trspos times the diagonal entry, the corresponding node is
!         transferred unaggregated to the coarse grid
!  kaptg_ is the threshold used to accept or not a tentative aggregate
!         when applying the coarsening algorithms from [2,3];
!         kaptg_blocdia is used for control based on bloc diagonal smoother [2];
!         kaptg_dampJac is used for control based on Jacobi smoother [3].
!  checkdd is the threshold to keep outside aggregation nodes where
!         the matrix is strongly diagonally dominant (based on mean of row
!         and column);
!         In fact, AGMG use the maximum of |checkdd| and of the default value
!            according to kaptg_ as indicated in [2,3]
!            (hence |checkdd| < 1 ensures that one uses this default value)
!         checkdd <0 : consider |checkdd|, but base the test on minus the
!                sum of offdiagonal elements, without taking the absolute value.
!  targetcoarsefac is the target coarsening factor (parameter tau in the main
!         coarsening algorithm in [2,3]): futher pairwise aggregation passes
!         are omitted once the number of nunzero entries has been reduced by a
!         factor of at least targetcoarsefac.
!  fracnegrcsum: if, at some level, more than fracnegrcsum*n nodes have
!         negative mean row and column sum, then the aggregation algorithm
!         of [2,3] is modified, exchanging all diagonal entries for the mean
!         row and column sum (that is, the algorithm is applied to a
!         modified matrix with mean row and colum sum enforced to be zero)
    REAL(kind(0.0d0)) :: trspos=0.45
    REAL(kind(0.0d0)) :: kaptg_blocdia=8, kaptg_dampJac=10
    REAL(kind(0.0d0)) :: checkdd=-0.5
    REAL(kind(0.0d0)) :: targetcoarsefac=4
    REAL(kind(0.0d0)) :: fracnegrcsum=0.25
!!!!!!!!!!!!!!!!!!!! END of PARAMETERS DEFINITION -----------------
!!!!!!!!!!!!!!!!!!! Internal variables declaration
!
    TYPE InnerData
       REAL(kind(0.0d0)), DIMENSION(:), POINTER :: a
       INTEGER, DIMENSION(:), POINTER :: ja
       INTEGER, DIMENSION(:), POINTER :: ia
       INTEGER, DIMENSION(:), POINTER :: il
       INTEGER, DIMENSION(:), POINTER :: iu
       REAL(kind(0.0d0)), DIMENSION(:), POINTER :: p
       INTEGER, DIMENSION(:), POINTER :: idiag
       INTEGER, DIMENSION(:), POINTER :: ind
       INTEGER, DIMENSION(:), POINTER :: iext
       INTEGER, DIMENSION(:), POINTER :: ilstout
       INTEGER, DIMENSION(:), POINTER :: lstout
       INTEGER, DIMENSION(:), POINTER :: ilstin
       INTEGER, DIMENSION(:), POINTER :: lstin
       INTEGER, DIMENSION(:), POINTER :: iblockl
    END TYPE InnerData
!
    TYPE(InnerData) :: dt(maxlev)
    REAL(kind(0.0d0)), ALLOCATABLE :: scald(:)
    INTEGER :: nn(maxlev),kstat(2,maxlev)=0,innermax(maxlev)
    INTEGER :: nlev,nwrkcum,iout,nrst,nbblock
    INTEGER :: maxcoarset,maxcoarseslowt
    REAL(kind(0.0d0)) :: memi=0.0,memax=0.0,memr=0.0,mritr,rlenilen,smoothtp,omeg
    REAL(kind(0.0d0)) :: wcplex(4),fracnz(maxlev),flop=0.0d0
    LOGICAL :: spd,wfo,wff,woo,allzero,trans,transint,zerors,gcrin
    LOGICAL :: coasing,coarsit
    REAL(kind(0.0d0)), PARAMETER :: cplxmax=3.0, xsi=0.6
    REAL(kind(0.0d0)), PARAMETER :: repsmach=SQRT(EPSILON(1.0d0)),epsmach=EPSILON(1.0d0)
    REAL(kind(0.0d0)), PARAMETER :: UNP=1.0d0+10*epsmach
    INTEGER :: nlc(2),nlcp(2),nlc1(2),icum,npassr
    INTEGER :: imult
    REAL(kind(0.0d0)) :: ngl(2),nglp(2),nlctot(2),ngl1(2),ngltot(2),RELRESL1
    INTEGER, DIMENSION(:), POINTER :: iextL1,ilstoutL1,lstoutL1,ilstinL1
    INTEGER, DIMENSION(:), POINTER :: lstinL1,ipermL1
    INTEGER, PARAMETER :: IRANK=-9999
    REAL(kind(0.0d0)) :: checkddJ, checkddB
  END MODULE dag2l_mem
!-----------------------------------------------------------------------
! END of Internal variables declaration
!-----------------------------------------------------------------------
MODULE dag2l_AGGROUTINES
CONTAINS
  SUBROUTINE dag2l_aggregation(l,n,a,ja,ia,nc,listrank)
!
!
!
!
!
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: l,n,nc
    INTEGER :: ja(*),ia(n+1)
    INTEGER, OPTIONAL :: listrank(*)
    REAL (kind(0.0d0)) :: a(*), dum
    INTEGER :: ier,i,j,k,maxdg,np,kpass,nzc,m1,ndd,nzp,isize,nddp,npass1,nz,i0
    LOGICAL :: skipass
    !
    INTEGER, POINTER, DIMENSION(:) :: jan,ian,idiagn,iextn,ind2,lcg,lcgn
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: an
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sinn
    !
    INTEGER, POINTER, DIMENSION(:) :: jap,iap,idiagp,iextp,lcg1
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sip
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ldd,iw,iperm,riperm
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: si1,w
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: wc
    !
    CHARACTER (len=80) :: lout
    !
    IF (l .EQ. 1) THEN
       IF (wfo) THEN
          WRITE (lout,901) IRANK
          call mexPrintf(lout//char(10))
       END IF
       IF (wff) THEN
          IF (spd) THEN
             WRITE (lout,906)
             call mexPrintf(lout//char(10))
          ELSE IF ( .NOT. transint) THEN
             WRITE (lout,908)
             call mexPrintf(lout//char(10))
          END IF
          IF (.not.spd) then
             WRITE (lout,902) 'Jacobi',kaptg_dampJac,checkddJ
             call mexPrintf(lout//char(10))
          ELSE
             WRITE (lout,902) 'BlockD',kaptg_blocdia,checkddB
             call mexPrintf(lout//char(10))
          END IF
          IF (checkdd < 0) THEN
             WRITE (lout,904)
             call mexPrintf(lout//char(10))
          END IF
          WRITE(lout,903) npass,targetcoarsefac
          call mexPrintf(lout//char(10))
          WRITE (lout,905) trspos
          call mexPrintf(lout//char(10))
       END IF
    END IF
    !
    IF (l .EQ. 1) THEN
       ALLOCATE(si1(n),ind2(n),iperm(n),riperm(n))
       memi=memi+4*n+1
       memr=memr+n
       CALL dag2l_setCMK(n,ja,ia,dt(l)%idiag,riperm,iperm)
    ELSE
       ALLOCATE(si1(n),ind2(n),iperm(n))
       memi=memi+2*n
       memr=memr+n
       iperm(1:n)=1
    END IF
    memax=MAX(memax,memr+memi*rlenilen)
    !
    call dag2l_prepareagg(n,a,ja,ia,dt(l)%idiag,ind2,iperm,si1,ndd,l )
    !
       IF (wfo) THEN
          IF (zerors) THEN
             WRITE (lout,907) IRANK
             call mexPrintf(lout//char(10))
          END IF
907       FORMAT(i3,'*  Too much nodes with neg. row+col. sum:',  &
           ' diag. modified for aggregation')
       END IF
    IF (ndd .EQ. n) THEN
       nc=0
       nzc=0
       DEALLOCATE(si1,iperm,ind2)
       memi=memi-2*n
       memr=memr-n
       IF (l .EQ. 1) THEN
          DEALLOCATE(riperm)
          memi=memi-n
          IF (smoothtype .EQ. 0) omeg=1.0d0
       END IF
       GOTO 999
    END IF
    !
    ALLOCATE(ldd(ndd),lcg(2*(n-ndd)))
    memi=memi+2*n-ndd
    memax=MAX(memax,memr+memi*rlenilen)
    !
    IF (dble(n) .GT. targetcoarsefac*(n-ndd)) THEN
       skipass=.TRUE.
       npass1=npass+1
    ELSE
       skipass=.FALSE.
       npass1=npass
    END IF
    !
    !
    IF (l > 1) THEN
       IF (spd) THEN
          CALL dag2l_findpairs_SI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       ELSE
          CALL dag2l_findpairs_GI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       END IF
       DEALLOCATE(iperm)
       memi=memi-n
    ELSE
       IF (spd) THEN
          CALL dag2l_findpairs_SI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       ELSE
          CALL dag2l_findpairs_GI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       END IF
       DEALLOCATE(iperm,riperm)
       memi=memi-2*n
    END IF
10  CONTINUE
    nz=ia(n+1)-1
    !
    IF (npass1.GT.1) THEN
       isize=nc
    ELSE
       isize=1
    END IF
    ALLOCATE(an(nz-2*(n-nc)+ndd),jan(nz-2*(n-nc)+ndd)             &
         ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
    memi=memi+nz-2*(n-nc)+ndd+2*(nc+1)+2*nc+isize
    memr=memr+nz-2*(n-nc)+ndd+nc
    memax=MAX(memax,memr+memi*rlenilen)
    CALL dag2l_setcg(n,a,ja,ia,dt(l)%idiag,si1,ind2,lcg,nc,an,jan,ian      &
         ,idiagn,sinn,npass1.GT.1,maxdg,iw,wc )
    DEALLOCATE(wc,iw)
    memi=memi-2*nc
    memr=memr-nc
    nzc=ian(nc+1)-1
    IF (dble(nz).GT.targetcoarsefac*nzc .OR. npass1.LE.1) THEN
       DEALLOCATE(si1,sinn,lcg)
       IF(ALLOCATED(ldd)) DEALLOCATE(ldd)
       memi=memi-2*(n-ndd)
       memr=memr-n-isize
       dt(l)%ind  => ind2
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(ind2,an,jan,ian,idiagn)
       GOTO 999
    END IF
    !
    DEALLOCATE(ind2)
    memi=memi-n
    !
    lcg1 => lcg
    NULLIFY(lcg)
    m1=1
    !
    !
    !
    DO kpass=2,npass1
       m1=2*m1
       np  = nc
       nzp = nzc
       ap     => an
       jap    => jan
       iap    => ian
       idiagp => idiagn
       sip    => sinn
       NULLIFY(an,jan,ian,idiagn,sinn)
       ALLOCATE(lcg(2*np),ind2(np),w(maxdg),iw(maxdg))
       memi=memi+maxdg+3*np
       memr=memr+maxdg
       memax=MAX(memax,memr+memi*rlenilen)
       !
       !
       ind2(1:np)=-1
       IF (spd) THEN
          CALL dag2l_findpairs_SF(l,np,ap,jap,iap,idiagp,sip  &
               ,ind2,lcg,nc                                  &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       ELSE
          CALL dag2l_findpairs_GF(l,np,ap,jap,iap,idiagp,sip     &
               ,ind2,lcg,nc                                     &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       END IF
       DEALLOCATE(w,iw)
       memi=memi-maxdg
       memr=memr-maxdg
       IF (kpass.NE.npass1) THEN
          isize=nc
       ELSE
          isize=1
          DEALLOCATE(si1)
          memr=memr-n
       END IF
       !
       !
       ALLOCATE(an(nzp-2*(np-nc)),jan(nzp-2*(np-nc))                 &
            ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
       memi=memi+nzp-2*(np-nc)+2*(nc+1)+2*nc+isize
       memr=memr+nzp-2*(np-nc)+nc
       memax=MAX(memax,memr+memi*rlenilen)
       !
       !
       CALL dag2l_setcg(np,ap,jap,iap,idiagp,sip,ind2,lcg,nc,an     &
            ,jan,ian,idiagn,sinn,kpass.NE.npass1,maxdg,iw,wc )
       memi=memi-SIZE(jap)-2*(np+1)-np-2*nc
       memr=memr-SIZE(ap)-SIZE(sip)-nc
       DEALLOCATE(ap,jap,iap,idiagp,sip,ind2,wc,iw)
       !
       !
       ALLOCATE(lcgn(2*m1*nc))
       memi=memi+2*m1*nc
       memax=MAX(memax,memr+memi*rlenilen)
       CALL dag2l_lcgmix(nc,m1,lcg1,lcg,lcgn)
       memi=memi-SIZE(lcg)-SIZE(lcg1)
       DEALLOCATE(lcg,lcg1)
       lcg1 => lcgn
       NULLIFY(lcgn)
       nzc=ian(nc+1)-1
       IF ( kpass.NE.npass1 .AND. dble(nz).GT.targetcoarsefac*nzc ) THEN
          DEALLOCATE(si1)
          memr=memr-n
          EXIT
       END IF
    END DO
    !
    memr=memr-SIZE(sinn)
    DEALLOCATE(sinn)
    !
    ALLOCATE(dt(l)%ind(n))
    memi=memi+n
    memax=MAX(memax,memr+memi*rlenilen)
    CALL dag2l_setind(nc,ndd,ldd,lcg1,2*m1,dt(l)%ind)
    memi=memi-ndd-SIZE(lcg1)
    DEALLOCATE(lcg1,ldd)
    !
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(an,jan,ian,idiagn)
999 CONTINUE
    !
    RETURN
901 FORMAT(i3,'*SETUP: Coarsening by multiple pairwise aggregation')
902 FORMAT('****       Quality threshold (',A6,'):',f6.2, &
         ' ;  Strong diag. dom. trs:',f5.2)
903 FORMAT('****         Maximal number of passes:',i3,     &
         '  ; Target coarsening factor:',f5.2)
904 FORMAT('****           Diag. dom. checked w.r.t. sum of offdiag', &
          ' (no absolute vaues)')
905 FORMAT('****',22x,'Threshold for rows with large pos. offdiag.:',f5.2)
906 FORMAT('****  Rmk: Setup performed assuming the matrix symmetric')
908 FORMAT('****  Rmk: Setup performed for the transpose of the input matrix')
  END SUBROUTINE dag2l_aggregation
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_setCMK(n,ja,ia,idiag,riperm,iperm,iext)
    !
    !
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),riperm(*),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    LOGICAL :: exc
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,mindg,kdim
    INTEGER :: ifirst,ilast,kb,i0,i2rec,i2rec0
    CHARACTER (len=80) :: lout
    !
    ifirst=1
    ilast=n
    i2=ifirst
    i2rec=i2-1
    mindg=n+1
    DO i = ifirst,ilast
       dg=ia(i+1)-ia(i)
       IF (dg .GT. 1) THEN
          iperm(i)=-dg
          IF (dg.LT.mindg) THEN
             mindg=dg
             jj=i
          END IF
       ELSE
          riperm(i2)=i
          iperm(i)=i2
          i2=i2+1
       END IF
    ENDDO
    IF (i2 .GT. i2rec+1) THEN
       IF (wfo) THEN
          WRITE(lout,901) IRANK,i2-i2rec-1
       END IF
       i2rec=i2-1
    END IF
    i2rec0=i2rec
    !
    ijs=ifirst-1
    i1=i2
15  CONTINUE
    !
    IF (i2 .LE. ilast) THEN
      riperm(i2)=jj
      iperm(jj)=i2
    END IF
    !
    DO WHILE (i1.LE.i2 .AND. i2.LT.ilast)
       !
       !
       i=riperm(i1)
       ijs1=i2+1
       !
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       DO kk = j1,jd-1
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       DO kk = jd+1,j2
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       !
       ijs2=i2
       exc=.TRUE. .AND. ijs2.GT.ijs1
       DO WHILE(exc)
          exc=.FALSE.
          DO kk=ijs1+1,ijs2
             IF( iperm(riperm(kk)) .GT. iperm(riperm(kk-1)) )THEN
                j=riperm(kk)
                riperm(kk)=riperm(kk-1)
                riperm(kk-1)=j
                exc=.TRUE.
             END IF
          END DO
       END DO
       DO kk=ijs1,ijs2
          iperm(riperm(kk))=kk
       END DO
       !
       i1=i1+1
    END DO
    IF (i2 .LT. ilast .OR. i2rec .GE. ifirst) THEN
       IF (wfo) THEN
          IF (i2-i2rec > 20) THEN
            IF (i2rec > i2rec0) THEN
              WRITE(lout,903) IRANK,i2rec-i2rec0
              call mexPrintf(lout//char(10))
            ENDIF
            WRITE(lout,902) IRANK,i2-i2rec
            call mexPrintf(lout//char(10))
            i2rec0=i2
          ELSE IF (i2 .GE. ilast) THEN
            WRITE(lout,903) IRANK,i2-i2rec0
          END IF
          i2rec=i2
       END IF
    END IF
    IF (i2 .LT. ilast) THEN
       !
       jj=0
       DO WHILE (jj .EQ. 0)
          ijs=ijs+1
          IF (ijs .GT. ilast) THEN
             mindg=mindg+1
             ijs=ifirst
          END IF
          ijs1=ijs
          IF (iperm(ijs1).LT.0 .AND. ia(ijs1+1)-ia(ijs1).EQ.mindg) &
               jj=ijs1
       END DO
       i2=i2+1
       GOTO 15
    END IF
    !
 !
 !
    RETURN
901 FORMAT(i3,'*',i12,' row(s) without offdiagonal connection')
902 FORMAT(i3,'*',i12,' unknowns in this connected component')
903 FORMAT(i3,'*',i12,' unknowns in small connected component(s) (size<=20)')
  END SUBROUTINE dag2l_setCMK
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_prepareagg(n,a,ja,ia,idiag,ind2,iperm,si,ndd,l,iext)
    !
    !
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind2(n),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER :: ndd, l
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) , TARGET :: si(n)
    REAL(kind(0.0d0)) :: checkddl,oda,odm,ods,vald
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,kdim,nnegrcs
    INTEGER :: ifirst,ilast,i0,kb
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: odmax,odabs,osi
    !
    IF (.NOT.spd) THEN
       checkddl=checkddJ
    ELSE
       checkddl=checkddB
    END IF
    !
    IF (.NOT.spd) THEN
       !
       !
       IF (checkdd > 0) THEN
          ALLOCATE(odmax(n),odabs(n))
          memi=memi+2*n
          odabs(1:n)=0.0d0
       ELSE
          ALLOCATE(odmax(n))
          memi=memi+n
       END IF
       osi => si(1:n)
       si(1:n)=0.0d0
       odmax(1:n)=0.0d0
       memax=MAX(memax,memr+memi*rlenilen)
       ifirst=1
       ilast=n
        DO i=ilast,ifirst,-1
          j =ia(i)
          jd=idiag(i)
          jj=ia(i+1)-1
          DO k = j,jd-1
            jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
          DO k = jd+1,jj
            jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
        ENDDO
    END IF
    !
    ndd=0
    nnegrcs=0
    !
    ifirst=1
    ilast=n
     DO i=ifirst,ilast
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       vald = a(jd)
       odm=0.0d0
       oda=0.0d0
       ods=0.0d0
       DO kk = j1,jd-1
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       DO kk = jd+1,j2
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       !
       IF (.NOT.spd) THEN
          ods=(osi(i)+ods)/2
          odm=max(odm,odmax(i))
          IF (checkdd > 0) oda=(oda+odabs(i))/2
       END IF
       !
       IF ((vald+ods) .LT. -repsmach*ABS(vald)) nnegrcs=nnegrcs+1
       !
       !
       si(i)=-ods
       IF ( (checkdd.GT.0 .AND. vald.GT.checkddl*oda)     &
            .OR. (checkdd.LT.0 .AND. vald.GT.checkddl*abs(ods)) ) THEN
          !
          ind2(i)=0
          ndd=ndd+1
       ELSE
          ind2(i)=-1
          IF (odm .GT. trspos*vald) iperm(i)=0
       ENDIF
     END DO
    !
    !
    zerors=.FALSE.
    IF (nnegrcs.GT.fracnegrcsum*n) THEN
       zerors=.TRUE.
       ndd=0
       ind2(1:n)=-1
    END IF
    !
    IF (.NOT.spd) THEN
       IF (checkdd > 0) THEN
          DEALLOCATE(odmax,odabs)
          memi=memi-2*n
       ELSE
          DEALLOCATE(odmax)
          memi=memi-n
       END IF
    END IF
    RETURN
  END SUBROUTINE dag2l_prepareagg
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_findpairs_GF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
!
!
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    REAL(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    CHARACTER (len=80) :: lout
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       CALL dag2l_checktentagg_GF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_checktentagg_GF
!
!
!
!
!
    INTEGER (KIND=selected_int_kind(12)), PARAMETER :: mm=max(2**(10),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m, setdim, l2, k2
    INTEGER (KIND=selected_int_kind(12)) :: info, setdim1
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)/2
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
       DO k=ia1(jj),idiag1(jj)-1
          DO l=1,j-1
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)/2
                W(j,l)=W(j,l)+alpha
                W(l,j)=W(j,l)
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
       v(j)=W(j,j)
   !
   !
   !
       W(j,j)=umdbndmum1*W(j,j)-abs(sig(j))
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
!
    beta=dbndmum1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
!
    acc=.FALSE.
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE dag2l_checktentagg_GF
  END SUBROUTINE dag2l_findpairs_GF
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_findpairs_SF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
!
!
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    REAL(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    CHARACTER (len=80) :: lout
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          epsr=repsmach*ABS(vals)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < -epsr) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < -epsr) CYCLE
          IF (vals > 0.0d0) THEN
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       CALL dag2l_checktentagg_SF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_checktentagg_SF
!
!
!
!
!
    INTEGER (KIND=selected_int_kind(12)), PARAMETER :: mm=max(2**(10),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m, setdim, l2, k2
    INTEGER (KIND=selected_int_kind(12)) :: info, setdim1
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
   !
       W(j,j)=W(j,j)-abs(sig(j))
   !
   !
   !
   !
       tmp=2*abs(sig(j))
       W(j,j)=W(j,j)-bndmum1m1*tmp
   !
       v(j)=tmp+AGe(j)
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
!
    beta=bndmum1m1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
!
    acc=.FALSE.
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE dag2l_checktentagg_SF
  END SUBROUTINE dag2l_findpairs_SF
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_findpairs_GI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
!
!
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    CHARACTER (len=80) :: lout
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  END SUBROUTINE dag2l_findpairs_GI
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_findpairs_SI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
!
!
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    CHARACTER (len=80) :: lout
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          epsr=repsmach*ABS(vals)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < -epsr) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < -epsr) CYCLE
          IF (vals > 0.0d0) THEN
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  END SUBROUTINE dag2l_findpairs_SI
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_findpairs_GI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
!
!
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    CHARACTER (len=80) :: lout
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  END SUBROUTINE dag2l_findpairs_GI1
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_findpairs_SI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
!
!
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    CHARACTER (len=80) :: lout
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          epsr=repsmach*ABS(vals)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < -epsr) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < -epsr) CYCLE
          IF (vals > 0.0d0) THEN
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  END SUBROUTINE dag2l_findpairs_SI1
!------------------------------------------------------------------
  SUBROUTINE dag2l_lcgmix(nc,m,lcg1,lcg,lcgn)
    INTEGER :: nc,m,lcg1(m,*),lcg(2,*),lcgn(2*m,*),i,l,l1,l2
    IF (m.eq.2) THEN
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(4,i)=-1
          ELSE IF(lcg(2,i) .LT. 0) THEN
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=-1
                lcgn(4,i)=-1
             ELSE
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=lcg1(2,lcg(1,i))
                lcgn(4,i)=-2
             END IF
          ELSE
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-2
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(3,i)=lcg1(2,lcg(2,i))
                   lcgn(4,i)=-3
                END IF
             ELSE
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-3
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=lcg1(2,lcg(2,i))
                END IF
             END IF
          END IF
       END DO
    ELSE
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(2*m,i)=-1
          ELSE
             lcgn(2,i)=-1
             l1=m
             IF (lcg1(m,lcg(1,i)).LT.0) l1=-lcg1(m,lcg(1,i))
             lcgn(1:l1,i)=lcg1(1:l1,lcg(1,i))
             IF(lcg(2,i) .LT. 0) THEN
                l=l1
             ELSE
                l2=m
                IF (lcg1(m,lcg(2,i)).LT.0) l2=-lcg1(m,lcg(2,i))
                lcgn(l1+1:l1+l2,i)=lcg1(1:l2,lcg(2,i))
                l=l1+l2
             END IF
             IF(l .LT. 2*m) lcgn(2*m,i)=-l
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE dag2l_lcgmix
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_setind(nc,ndd,ldd,lcg,m,ind)
    INTEGER :: nc,m,lcg(m,*),nll,ldd(ndd),ind(*),i,k,l
    DO i=1,ndd
       ind(ldd(i))=0
    END DO
    DO i=1,nc
       l=m
       IF (lcg(m,i) .LT. 0) l=-lcg(m,i)
       DO k=1,l
          ind(lcg(k,i))=i
       END DO
    END DO
    RETURN
  END SUBROUTINE dag2l_setind
!-----------------------------------------------------------------------
  SUBROUTINE dag2l_setcg(n,a,ja,ia,idiag,si,ind,lcg           &
       ,nc,a2,ja2,ia2,idiag2,si2,ysi,maxdg,iw,w,iext,iext2)
    USE dag2l_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),nc,lcg(2,nc)
    INTEGER :: ja2(*),ia2(nc+1),idiag2(nc),maxdg
    INTEGER, TARGET :: iw(2*nc)
    INTEGER, OPTIONAL :: iext(*),iext2(*)
    REAL(kind(0.0d0)) :: a(*),a2(*),w(nc),vald
    REAL(kind(0.0d0)) :: si(n),si2(*),sii
    LOGICAL :: ysi
    INTEGER :: nz,nzu,i,jj,jc,jcol,ki,kb,kf,jpos
    INTEGER, POINTER, DIMENSION(:) :: iw1, iw2
    !
    iw1 => iw(1:nc)
    iw2 => iw(nc+1:2*nc)
    !
    nz = 0
    iw1(1:nc)=0
    maxdg=0
    ia2(1)=1
    DO i = 1,nc
       sii=0.0d0
       vald=0.0d0
       nzu=0
       DO ki= 1,2
          jj = lcg(ki,i)
          IF (ki.EQ.1 .OR. jj.GT.0) THEN
             IF (ysi) sii=sii+si(jj)
             kf = ia(jj+1) - 1
             DO kb = ia(jj),kf
                jc = ja(kb)
                jcol = ind(jc)
                IF (jcol .GT. 0) THEN
                   IF (jcol .LT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nz = nz+1
                         ja2(nz) = jcol
                         iw1(jcol) = nz
                         a2(nz) = a(kb)
                      ELSE
                         a2(jpos) = a2(jpos) + a(kb)
                      ENDIF
                   ELSE IF (jcol .GT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nzu = nzu+1
                         iw2(nzu) = jcol
                         iw1(jcol) = nzu
                         w(nzu) = a(kb)
                      ELSE
                         w(jpos) = w(jpos) + a(kb)
                      ENDIF
                   ELSE
                      vald=vald+a(kb)
                      IF (ysi .AND. jc.NE.jj) sii=sii+a(kb)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       nz=nz+1
       a2(nz)=vald
       idiag2(i)=nz
       ja2(nz)=i
       a2(nz+1:nz+nzu)=w(1:nzu)
       ja2(nz+1:nz+nzu)=iw2(1:nzu)
       nz=nz+nzu
       maxdg=max(maxdg,nz-ia2(i))
       DO kb = ia2(i), nz
          iw1(ja2(kb))=0
       ENDDO
       IF (ysi) si2(i)=sii
       ia2(i+1)=nz+1
    ENDDO
    RETURN
  END SUBROUTINE dag2l_setcg
!-----------------------------------------------------------------------
END MODULE dag2l_AGGROUTINES
SUBROUTINE dag2l_twolev(l,n,a,ja,ia,ind,sym,npas,kap,checkd,targetcf  &
                        ,fracnzrs,trspo,iprint )
  USE dag2l_mem
  USE dag2l_AGGROUTINES
  IMPLICIT NONE
  INTEGER :: l, l1, n, ja(*), ia(n+1), ind(n), sym, npas, iprint
  REAL(kind(0.0d0)) :: a(*), kap, checkd, targetcf, fracnzrs, trspo
  INTEGER :: nc, i, k
  transint=.TRUE.
  l1=1
  IF (l.NE.1) l1=2
  spd=sym.EQ.1
  npass=min(npas,10)
  kaptg_blocdia=kap
  kaptg_dampJac=kap
  targetcoarsefac=targetcf
  fracnegrcsum=fracnzrs
  trspos=trspo
  checkdd=checkd
  checkddJ=MAX(ABS(checkdd),kaptg_dampJac/(kaptg_dampJac-2))
  checkddB=MAX(ABS(checkdd),(kaptg_blocdia+1)/(kaptg_blocdia-1))
    wfo=.TRUE.
    iout=iprint
    IF (iprint <= 0) THEN
       iout=6
       wfo=.FALSE.
    ELSE IF (iprint == 5) THEN
       iout=6
    END IF
    wff=wfo
    !
    !
       ALLOCATE(dt(l1)%idiag(n+1))
       DO i=1,n
          DO k=ia(i),ia(i+1)-1
             IF (ja(k) .EQ. i) THEN
                dt(l1)%idiag(i)=k
                EXIT
             END IF
          END DO
       END DO
  imult=1
  CALL dag2l_aggregation(l1,n,a,ja,ia,nc)
  IF (nc > 0) THEN
    ind(1:n)=dt(l1)%ind(1:n)
    DEALLOCATE(dt(l1+1)%a,dt(l1+1)%ja,dt(l1+1)%ia,dt(l1+1)%idiag      &
              ,dt(l1)%ind)
  ELSE
    ind(1:n)=0
  END IF
  DEALLOCATE(dt(l1)%idiag)
  RETURN
END SUBROUTINE dag2l_twolev
!-----------------------------------------------------------------------
