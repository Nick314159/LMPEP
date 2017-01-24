!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                       MODULE EIGENSOLVE                             !!!
!!!!!              by D.A. Bini, L. Gemignani, F. Tisseur                 !!!
!!!!!                              v. 1.1                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Latest revision: March 2004  
! The program solves the generalized eigenvalue problem 
!                              det(T-xS)=0 
! by means of  Divide and Conquer  and Aberth's iteration where T is real
! symmetric tridiagonal and S is real diagonal. The problem must be 
! normalized so that the upper and lower diagonal entries of 
! T are 1 (normalized form). See the subroutine normalize.
!---------------------------------------------------------------------------
! THE ALGORITHM:
! The problem is split into two subproblems:
!  1)  det(T_1-x S_1)=0, where T_1 and S_1 are the leading
!   submatrices of size m=n/2 and
!  2)  det(T_2-x S_2)=0, where T_2 and S_2 are the trailing submatrices
!  of size n-m
! The original problem det(T-xS)=0 is solved by means of Aberth's
! iteration choosing as initial approximations the solution of the 
! two subproblems slightly perturbed so that there cannot be double 
! solutions.
! The two subproblems are recursively solved by means of the same technique.
! The computation of the ratio p(x)/p'(x) needed by Aberth's method, where
! p(x)=det(T-xS), is performed by means of QR factorization of the 
! matrix T-xS.
!----------------------------------------------------------------------------
! THE SUBROUTINES:
! EIGEN:    recursive subroutine 
! ABERTH:   Aberth's iteration
! NEWTCORR: computes the Newton correction
! VALIDATE: computes inclusion radii for all the approximations
! POLY: evaluates det(T-x S)
! NORMALIZE: transform the generalized eigenvalue problem T-xS into
!         normalized form
! Quadruple precision subroutines:
! QABERTH: It refines the approximations in quadruple precision
! QNEWTCORR:
! QNORMALIZE: 
!---------------------------------------------------------------------------
! PARAMETERS
! maxit: maximum number of Aberth's iterations in the last recursive stage
! maxit_int: maximum number of Aberth's iterations in the intermediate stages
! eps: Machine epsilon
! correct: max value of the perturbation introduced in the merged eigenvalues
! debug: If set to .TRUE. then debug information is output in the file 
!        fort.11
!---------------------------------------------------------------------------
! CALLING LINE:
! CALL eigen(n,a,s,z,cond)
!
! n: size of the matrix
! a: diagonal elements of T
! s: diagonal elements of S
! z: complex vector with the approximations to the eigenvalues
! cond: real vector with the condition number estimates of the eigenvalues
!---------------------------------------------------------------------------

MODULE eigensolve
IMPLICIT NONE

INTEGER,  PARAMETER    :: dp = SELECTED_REAL_KIND(15, 60)
INTEGER,  PARAMETER    :: qp = SELECTED_REAL_KIND(30, 60)

INTEGER,  PARAMETER    :: maxit = 500,  maxit_int=20
REAL(DP), PARAMETER    :: eps=EPSILON(1.d0)
COMPLEX(DP), PARAMETER :: correct=(0,1.d0)*eps*10
LOGICAL, PARAMETER     :: debug=.false. 
PRIVATE                :: maxit, maxit_int, correct, eps, debug
! dp:  kind double precision
! maxit: maximum number of Aberth's iterations in the final step
! maxit_int: maximum number of Aberth's iters. in the intermediate stages
! correct: perturbation applied to the eigenvalues of the subproblems 
!        in the divide and conquer stage.

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                             SUBROUTINE VALIDATE                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preliminary version.
! Given the approximations z(i), i=1,..,n to the eigenvalues, the subroutine 
! computes the numbers
!        rad(i)=log_10(|det(T-z(i)S)|/|product_{j/=i}(z(i)-z(j)|
! such that the disks D_i, i=1,..,n, of center z(i) and radius 10**rad(i)
! are a set of Gerschgorin disks for the eigenvalues, i.e., their union
! contains all the eigenvalues and each connected component of this union, 
! made up by k disks, contains k eigenvalues.
! In this version the computation is not guaranteed, i.e., the cancelation
! in evaluating rad(i) may provide uncorrect bounds.
! The guaranteed version, relying on the running error analysis of the
! computation of det(T-xS), must be still written.
! INPUT:
!  n: size of the matrices T and S
!  a: real vector with the diagonal entries of the matrix T
!  s: real vector with diagonal entries of the matrix S
!  z: complex vector with the approximations of the eigenvalues
! OUTPUT:
!  rad: real vector with the logarithms to base 10 of the radii of the 
!  inclusion disks
  SUBROUTINE validate(n,a,s,z,rad)
    IMPLICIT NONE
    INTEGER                    :: n
    REAL(dp), DIMENSION(n)     :: a,s
    COMPLEX(dp), DIMENSION(n)  :: z
    REAL(dp),DIMENSION(n)      :: rad
    REAL(dp)::w
    !----------------------------------------------
    INTEGER                    :: i,j
    DO i=1,n
       CALL poly(n,a,s,z(i),rad(i))
       DO j=1,n
          IF(.NOT. j == i) THEN
             IF(ABS(z(i)-z(j))==0)THEN
                rad(i)=1.d300
                EXIT
             ELSE
                w=LOG(ABS(z(i)-z(j)))
                rad(i)=rad(i)-w+(n-j+1)*abs(w)*eps
             END IF
          END IF
       END DO
       DO j=1,n
          w=LOG(ABS(s(j)))
          rad(i)=rad(i)-w+abs(w)*(n-j+1)*eps
       END DO
       rad(i)=rad(i)-LOG(ABS(z(i)))+LOG(n*1.d0)
       rad(i)=rad(i)/LOG(10.0)
    END DO
  END SUBROUTINE validate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                             SUBROUTINE POLY                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! evaluates log(det(T-z S))
! INPUT:
!  n: size of the matrices T and S
!  aa: real vector with the diagonal entries of the matrix T
!  ss: real vector with diagonal entries of the matrix S
!  z : complex number
! OUTPUT
!  r: log(det(T-zS))
  SUBROUTINE POLY(N,AA,SS,Z,p)
    IMPLICIT NONE
    INTEGER               :: n
    REAL(dp),DIMENSION(n) :: aa, ss
    REAL(dp)              :: p
    COMPLEX(dp)           :: z
    !------------------------
    INTEGER                   :: i
    COMPLEX(dp), DIMENSION(n) :: a, al, r
    REAL(dp), DIMENSION(n)    :: be
    REAL(dp)                  :: sq, ar
    
    a=aa-z*ss
    ! Compute Givens rotations alfa beta
    ! Compute the entries of R:  r, s, t
    DO i=1,n-1
       sq=SQRT(ABS(a(i))**2+1.d0)
       be(i)=1.d0/sq
       al(i)=a(i)*be(i)
       r(i)=cnjg(al(i))*a(i)+be(i)
       IF(i==1)THEN
          a(i+1)=be(i)-al(i)*a(i+1)
       ELSE
          a(i+1)=-be(i)*al(i-1)-al(i)*a(i+1)
       END IF
    END DO
    r(n)=a(n)
    p=0.d0
    DO i=1,n
       ar=ABS(r(i))
       IF(ar==0)THEN
          ar=-300
       ELSE
          ar=LOG(ar)
          ar=max(ar,LOG(eps)) 
          ar=ar +abs(ar)*eps*(n-i+1)
       END IF
       p=p+ar
    END DO
  END SUBROUTINE POLY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                             SUBROUTINE EIGEN                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
  SUBROUTINE eigen(n,a,s,z,cond)
    IMPLICIT NONE
    INTEGER                    :: n
    REAL(dp), DIMENSION(n)     :: a,s,cond
    COMPLEX(dp), DIMENSION(n)  :: z
    !----------------------------------------------
    REAL(dp)                   :: iter
    INTEGER                    :: itermx
    itermx=0
    CALL reigen(n,a,s,z,cond,maxit,iter,itermx)
    WRITE(*,*)"average number of iterations per eigenvalue ", iter/n
    WRITE(*,*)"max number of iterations per eigenvalue     ",itermx
  END SUBROUTINE eigen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                             SUBROUTINE REIGEN                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recursive subroutine for the evaluation of the eigenvalues of T-x S
  RECURSIVE SUBROUTINE reigen(n,a,s,z,cond,mxit,iter,itermx)
    IMPLICIT NONE
    INTEGER                    :: n,mxit,itermx
    REAL(dp), DIMENSION(n)     :: a,s,cond
    COMPLEX(dp), DIMENSION(n)  :: z
    REAL(dp)                   :: iter
    !----------------------------------------------
    INTEGER                    :: i, nm, itermx1
    REAL(dp), DIMENSION(:)     :: aa, ss,cond1,cond2
    COMPLEX(dp), DIMENSION(:)  :: zz1, zz2
    REAL(dp)                   :: ca, cb, cc, discr
    COMPLEX(dp)                :: sqrtdiscr,crct,cru
    ALLOCATABLE                ::  aa, ss, zz1, zz2,cond1,cond2
    REAL(dp)                   :: iter1, iter2, iter3, nrm1, nrm2
    REAL                       :: ru1,ru2,ru
! matrix of size n=1
    IF(n==1)THEN
       z(1)=a(1)/s(1)
       iter=0
       RETURN
    END IF
! matrix of size n=2
    IF(n == 2)THEN
       iter=0;itermx=0
       ca=s(1)*s(2)
       cb=s(1)*a(2)+s(2)*a(1)
       cc=a(1)*a(2)-1.d0
       IF(cc==0)THEN
          z(1)=cb/ca
          z(2)=0.d0
       ELSE
          discr=cb**2-4*ca*cc
          IF(discr>0)THEN
             discr=SQRT(discr)
             sqrtdiscr=discr
          ELSE
             discr=SQRT(-discr)
             sqrtdiscr=discr*(0,1.d0)
          END IF
          z(1)=(cb-sqrtdiscr)/(2*ca)
          z(2)=(cb+sqrtdiscr)/(2*ca)
       END IF
       RETURN
    ELSE                               
       ! matrix of size n>2: split into two subproblems
       nm=(n+1)/2
       ALLOCATE(aa(nm), ss(nm))
       ALLOCATE(zz1(nm), zz2(n-nm),cond1(nm),cond2(n-nm))
       aa=a(1:nm)
!    aa(nm)=aa(nm)-1 ! for rank 2 correction
       ss=s(1:nm)
       nrm1=maxval((2+abs(aa))/abs(ss))
       ! first recursive call
       CALL reigen(nm,aa,ss,zz1,cond1,maxit_int,iter1,itermx)
       DEALLOCATE(aa, ss)
       ALLOCATE(aa(n-nm), ss(n-nm))
       aa=a(nm+1:n)
!    aa(1)=aa(1)-1 ! for rank 2 correction
       ss=s(nm+1:n)
       nrm2=maxval((2+abs(aa))/abs(ss))

       ! second recursive call
       CALL reigen(n-nm,aa,ss,zz2,cond2,maxit_int,iter2,itermx)
       ! slightly perturb the eigenvalues of the two subproblems
       ! for avoiding multiple roots and merge them together
       ! the correction has a relative perturbation of modulus at most 
       ! eps and an absolute perturbation of modulus at most eps*normA 
       crct=correct
       ru1=0.1;ru2=-0.2;
       DO i=1,nm
          CALL RANDOM_NUMBER(ru1)
          call RANDOM_NUMBER(ru2)
          cru=cmplx(0.5-ru1,0.5-ru2)*crct
          z(i)=zz1(i)*(1-cru)-cru*nrm1
       END DO
       crct=correct
       DO i=1,n-nm
          CALL RANDOM_NUMBER(ru1)
          call RANDOM_NUMBER(ru2)
          cru=cmplx(0.5-ru1,0.5-ru2)*crct
          z(nm+i)=zz2(i)*(1+cru)+cru*nrm2 
       END DO
       ! refine the approximations by means of Aberth's iteration
       iter=0;itermx=0
       CALL aberth(n,a,s,mxit,z,cond,iter3,itermx1)
       itermx=itermx1 !max(itermx,itermx1)
       iter=iter3   !iter=iter1+iter2+iter3
    END IF
  END SUBROUTINE reigen
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                            SUBROUTINE ABERTH                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs the Aberth iterations for refining the 
! approximation z to the eigenvalue
  SUBROUTINE aberth(n,a,s,mxit,z,cond,witer,nit)
    IMPLICIT NONE
    INTEGER                     :: n,mxit,nit
    REAL(dp), DIMENSION(n)      :: a, s, cond
    COMPLEX(dp), DIMENSION(n)   :: z, diag
    REAL(dp)                    :: witer
    !---------------------------------------------
    INTEGER                     :: i, j, nzer, global
    COMPLEX(dp)                 :: corr, abcorr, newt, nwt1,nwt2, den
    LOGICAL,DIMENSION(n)        :: cont
    REAL(dp)                    :: rat
    !---------------------------------------------
    cont=.TRUE.;  nzer=0; nit=0; global=0
    DO WHILE(nzer<n .AND. nit<=mxit)
       nit=nit+1
       DO i=1,n
          IF(cont(i)) THEN
             global=global+1
             ! compute Aberth's correction
             abcorr=0
             DO j=1,n
                IF(.NOT. j == i) THEN
                   IF(ABS(z(i)-z(j))==0)THEN
                      WRITE(11,*)"WARNING: zi=zj"
                   ELSE
                      abcorr=abcorr+1.d0/(z(i)-z(j))
                   END IF
                END IF
             END DO
             ! compute the Newton correction
             CALL newtcorr(n,a,s,z(i),newt,cont(i),cond(i))
             IF(.NOT. cont(i))THEN
                nzer=nzer+1
              END IF
             IF(debug) THEN
                IF(.NOT. CONT(I))THEN
                   WRITE(11,*)"n=",n,"nzer=",nzer,"nit=",nit,"z=",z(i)
                END IF
                IF(nit>=maxit-3)THEN
                   WRITE(11,*)"n=",n,"nzer=",nzer,"nit=",nit
                   WRITE(11,*)"z=",z(i),"abscorr=",newt/(1.d0-newt*abcorr),"newtcorr=",newt
                 END IF
             END IF
             den=1.d0-newt*abcorr
             IF(den==0)THEN
                corr=newt
                WRITE(11,*)"WARNING INF CORRECTION for n=",n
             ELSE
                corr=newt/den
             END IF
             z(i)=z(i)-corr
          END IF
       END DO
    END DO
   IF(nit>=maxit)THEN
      WRITE(*,*)"WARNING: Exceeded max number of iterations n=",n
      WRITE(11,*)"WARNING: Exceeded max number of iterations n=",n
   END IF
   IF(debug) THEN
      WRITE(11,*)"n=",n,"average iterations per root", (1.0*global)/n
   END IF
   witer=global
 END SUBROUTINE aberth

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                          SUBROUTINE NEWTCORR                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the Newton correction nwtc=p(z)/p'(z)
! where p(z)=det(T-zS) by means of the formula 
! nwtc=-1/trace( (S**(-1)T-zI)**(-1) )
! the diagonal entries of (T-zS)**-1 are computed by means of the
! QR factorization of T-zS. For the latter computation, the semiseparable
! structure of (T-zS)**-1 is used.
! It is assumed that T-zS is in normal form, i.e., the off-diagonal entries
! of T are 1.
! INPUT:
!  n  : size of the matrix T
!  aa : diagonal entries of T
!  ss : diagonal entries of S
!  z  : 
! OUTPUT
!  nwtc : Newton's correction
!  cont : .true. if the approximation z to the root of p is still poor
!         .false. if the approximation z is accurate
  SUBROUTINE newtcorr(n,aa,ss,z,nwtc,cont,rat)
    IMPLICIT NONE
    INTEGER               :: n
    REAL(dp),DIMENSION(n) :: aa, ss
    REAL(dp)              :: rat
    COMPLEX(dp)           :: z, nwtc
    LOGICAL               :: cont
    !------------------------
    INTEGER                   :: i
    COMPLEX(dp), DIMENSION(n) :: a, al, r, s, t, u, v, w
    REAL(dp), DIMENSION(n)    :: be
    REAL(dp)                  :: sq, mvr
    REAL(dp)                  :: mvl, anwtc,anw

    !------------------------    
    a=aa-z*ss
    mvl=MAXVAL(ABS(a))+2
    ! Compute Givens rotations alfa beta
    ! Compute the entries of R:  r, s, t
    DO i=1,n-1
       sq=SQRT(ABS(a(i))**2+1.d0)
       be(i)=1.d0/sq
       al(i)=a(i)*be(i)
       r(i)=cnjg(al(i))*a(i)+be(i)
       IF(i==1)THEN
          s(i)=cnjg(al(i))+be(i)*a(i+1)
          a(i+1)=be(i)-al(i)*a(i+1)
       ELSE
          s(i)=-cnjg(al(i))*al(i-1)+be(i)*a(i+1)
          a(i+1)=-be(i)*al(i-1)-al(i)*a(i+1)
       END IF
    END DO
    r(n)=a(n)
    mvr=MINVAL(ABS(r))
    IF(mvr<1.d-308)THEN
       cont=.FALSE.
       WRITE(11,*)"mvr=0", " n=",n, "z=",z
       RETURN
    END IF

    ! Scale R and compute scaled u, v
    s(1:n-1)=s(1:n-1)*be(1:n-1)
    t(1:n-2)=be(1:n-2)**2 *be(2:n-1)
    u(1)=1;
    u(2:n)=-al(1:n-1)
    v(n)=1
    DO i=1,n-1
       v(i)=cnjg(al(i))
    END DO
    
    ! solve the linear system
    w(n)=v(n)/r(n)
    w(n-1)=(v(n-1)-w(n)*s(n-1))/r(n-1)
    DO i=n-2,1,-1
       w(i)=(v(i)-w(i+1)*s(i)-w(i+2)*t(i))/r(i)
    END DO

    ! compute the trace
    w=w*u*ss
    nwtc=SUM(w)
    IF(nwtc==0)THEN
       WRITE(11,*)"WARNING : null trace n=",n," z=",z
       nwtc= 0
       RETURN
    END IF
    anw=abs(nwtc)
    rat=SUM(ABS(w))/anw  
    nwtc=-1.d0/nwtc
    anwtc=1.d0/anw
    IF(anwtc<mvl*rat*eps .or. anwtc<abs(z)*eps .or. anwtc<mvl*eps)THEN 
       cont=.FALSE.
    END IF
  END SUBROUTINE newtcorr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                         SUBROUTINE NORMALIZE                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine transforms the generalized tridiagonal eigenvalue problem 
! in normal form where T has off-diagonal entries equal to 1.
! INPUT:
!  n : size of the matrix
!  a : diagonal entries of T
!  b : upper diagonal entries of T
!  c : lower diagonal entries of T
!  s : diagonal entries of S
! OUTPUT:
!  a : diagonal entries of the normal form of T
!  s : diagonal entries of S in the normal form
  SUBROUTINE normalize(n,a,b,c,s)
    IMPLICIT NONE
    INTEGER :: n, i
    REAL(DP),DIMENSION(:)::a,b,c,s
    REAL(DP) :: t
    DO i=2,n-1
       t=1.d0/(c(i-1)*b(i-1))
       a(i)=a(i)*t
       s(i)=s(i)*t
       b(i)=b(i)/c(i-1)
       c(i)=c(i)/b(i-1)
    END DO
    t=1.d0/(c(n-1)*b(n-1))
    a(n)=a(n)*t
    s(n)=s(n)*t
  END SUBROUTINE normalize

  COMPLEX(dp) FUNCTION cnjg(x)
    IMPLICIT NONE
    COMPLEX(dp)::x
    cnjg=REAL(x)-AIMAG(x)*(0.d0,1.d0)
  END FUNCTION cnjg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                    QUADRUPLE PRECISION ROUTINES                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                            SUBROUTINE QABERTH                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs the Aberth iterations for refining the 
! approximation z to the eigenvalue in quadruple precision
  SUBROUTINE qaberth(n,a,s,z)
    IMPLICIT NONE
    INTEGER                     :: n,mxit,nit
    REAL(qp), DIMENSION(n)      :: a, s
    COMPLEX(qp), DIMENSION(n)   :: z, diag
    !---------------------------------------------
    INTEGER                     :: i, j, nzer, global
    COMPLEX(qp)                 :: corr, abcorr, newt, nwt1,nwt2, den
    LOGICAL,DIMENSION(n)        :: cont
    REAL(qp)                    :: rat,cond
    !---------------------------------------------
    cont=.TRUE.;  nzer=0; nit=0; global=0
    mxit=maxit
    DO WHILE(nzer<n .AND. nit<=mxit)
       nit=nit+1
       DO i=1,n
          IF(cont(i)) THEN
             global=global+1
             ! compute Aberth's correction
             abcorr=0
             DO j=1,n
                IF(.NOT. j == i) THEN
                   IF(ABS(z(i)-z(j))==0)THEN
                      WRITE(11,*)"WARNING: zi=zj"
                   ELSE
                      abcorr=abcorr+1.d0/(z(i)-z(j))
                   END IF
                END IF
             END DO
             ! compute the Newton correction
             CALL qnewtcorr(n,a,s,z(i),newt,cont(i),cond)
             IF(.NOT. cont(i))THEN
                nzer=nzer+1
              END IF
             IF(debug) THEN
                IF(.NOT. CONT(I))THEN
                   WRITE(11,*)"n=",n,"nzer=",nzer,"nit=",nit,"z=",z(i)
                END IF
                IF(nit>=maxit-3)THEN
                   WRITE(11,*)"n=",n,"nzer=",nzer,"nit=",nit
                   WRITE(11,*)"z=",z(i),"abscorr=",newt/(1.d0-newt*abcorr),"newtcorr=",newt
                 END IF
             END IF
             den=1.d0-newt*abcorr
             IF(den==0)THEN
                corr=newt
                WRITE(11,*)"WARNING INF CORRECTION for n=",n
             ELSE
                corr=newt/den
             END IF
             z(i)=z(i)-corr
          END IF
       END DO
    END DO
   IF(nit>=maxit)THEN
      WRITE(*,*)"WARNING: Exceeded max number of iterations n=",n
      WRITE(11,*)"WARNING: Exceeded max number of iterations n=",n
   END IF
   IF(debug) THEN
      WRITE(11,*)"n=",n,"average iterations per root", (1.0*global)/n
   END IF
   write(*,*)"max number of iterations per eigenvalue ",nit
 END SUBROUTINE qaberth
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                          SUBROUTINE QNEWTCORR                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the Newton correction nwtc=p(z)/p'(z)
! in quartuple precision
  SUBROUTINE qnewtcorr(n,aa,ss,z,nwtc,cont,rat)
    IMPLICIT NONE
    INTEGER               :: n
    REAL(qp),DIMENSION(n) :: aa, ss
    REAL(qp)              :: rat
    COMPLEX(qp)           :: z, nwtc
    LOGICAL               :: cont
    !------------------------
    INTEGER                   :: i
    COMPLEX(qp), DIMENSION(n) :: a, al, r, s, t, u, v, w
    REAL(qp), DIMENSION(n)    :: be
    REAL(qp)                  :: sq, mvr, qeps
    REAL(qp)                  :: mvl, anwtc,anw

    !------------------------    
    a=aa-z*ss
    mvl=MAXVAL(ABS(a))+2
    qeps=EPSILON(mvl)
    ! Compute Givens rotations alfa beta
    ! Compute the entries of R:  r, s, t
    DO i=1,n-1
       sq=SQRT(ABS(a(i))**2+1.d0)
       be(i)=1.d0/sq
       al(i)=a(i)*be(i)
       r(i)=qcnjg(al(i))*a(i)+be(i)
       IF(i==1)THEN
          s(i)=qcnjg(al(i))+be(i)*a(i+1)
          a(i+1)=be(i)-al(i)*a(i+1)
       ELSE
          s(i)=-qcnjg(al(i))*al(i-1)+be(i)*a(i+1)
          a(i+1)=-be(i)*al(i-1)-al(i)*a(i+1)
       END IF
    END DO
    r(n)=a(n)
    mvr=MINVAL(ABS(r))
    IF(mvr<1.d-308)THEN
       cont=.FALSE.
       WRITE(11,*)"mvr=0", " n=",n, "z=",z
       RETURN
    END IF

    ! Scale R and compute scaled u, v
    s(1:n-1)=s(1:n-1)*be(1:n-1)
    t(1:n-2)=be(1:n-2)**2 *be(2:n-1)
    u(1)=1;
    u(2:n)=-al(1:n-1)
    v(n)=1
    DO i=1,n-1
       v(i)=qcnjg(al(i))
    END DO
    
    ! solve the linear system
    w(n)=v(n)/r(n)
    w(n-1)=(v(n-1)-w(n)*s(n-1))/r(n-1)
    DO i=n-2,1,-1
       w(i)=(v(i)-w(i+1)*s(i)-w(i+2)*t(i))/r(i)
    END DO

    ! compute the trace
    w=w*u*ss
    nwtc=SUM(w)
    IF(nwtc==0)THEN
       WRITE(*,*)"WARNING : null trace n=",n," z=",z
       WRITE(11,*)"WARNING : null trace n=",n," z=",z
       nwtc= 0
       RETURN
    END IF
    anw=abs(nwtc)
    rat=SUM(ABS(w))/anw  
    nwtc=-1.d0/nwtc
    anwtc=1.d0/anw
    IF(anwtc<mvl*rat*qeps .or. anwtc<abs(z)*qeps .or. anwtc<mvl*qeps)THEN 
       cont=.FALSE.
    END IF
  END SUBROUTINE qnewtcorr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                         SUBROUTINE QNORMALIZE                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine transforms the generalized tridiagonal eigenvalue problem 
! in normal form where T has off-diagonal entries equal to 1.
! Computations are performed in quartuple precision
  SUBROUTINE qnormalize(n,a,b,c,s)
    IMPLICIT NONE
    INTEGER :: n, i
    REAL(QP),DIMENSION(:)::a,b,c,s
    REAL(QP) :: t
    DO i=2,n-1
       t=1.d0/(c(i-1)*b(i-1))
       a(i)=a(i)*t
       s(i)=s(i)*t
       b(i)=b(i)/c(i-1)
       c(i)=c(i)/b(i-1)
    END DO
    t=1.d0/(c(n-1)*b(n-1))
    a(n)=a(n)*t
    s(n)=s(n)*t
  END SUBROUTINE qnormalize


  COMPLEX(qp) FUNCTION qcnjg(x)
    IMPLICIT NONE
    COMPLEX(qp)::x
    qcnjg=REAL(x)-AIMAG(x)*(0.d0,1.d0)
  END FUNCTION qcnjg


END MODULE eigensolve




