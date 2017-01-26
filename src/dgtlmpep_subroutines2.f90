MODULE dgtlmpep_subroutines2
USE dslmpep_subroutines
IMPLICIT NONE
REAL(dp), PARAMETER :: two=2.0_dp
COMPLEX(dp), PARAMETER :: ctwo=DCMPLX(two)
CONTAINS

!************************************************************************
!			SUBROUTINE DGTLM				*
!************************************************************************
! Compute the eigenvalues and eigenvectors of the matrix polynomial p 	*
! with real general-tridiagonal coefficients of degree d and size n 	*
! using Laguerre's method. The norm of the coefficients of p are stored *
! in ncoeff, iseed is used when computing an upper bound on the backward* 
! error, stored in berr of each eigenpair. Eigenvalues are stored in 	*
! (er,ei) and eigenvectors in (xr,xi) and (yr,yi).			*
!************************************************************************
SUBROUTINE dgtlm(pdl, pd, pdu, er, ei, berr, ncoeff, iseed, d, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, n
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: ncoeff(*), pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: berr(*), er(*), ei(*)
!local scalars
LOGICAL :: check
INTEGER :: die, dze, i, it, td
REAL(dp) :: tol
!intrinsic procedures
INTRINSIC :: DABS, MAX

!initial estimates
CALL dgtstart(pdl, pd, pdu, er, ei, ncoeff, d, die, dze, n)
td=n*d-die
!Laguerre's method
DO i=dze+1,td
  DO it=1,itmax
    check=(it==itmax)
    tol=MAX(eps*DCMOD(er(i),ei(i)), eps)
    IF(DABS(ei(i))<tol) THEN
      !update real eigenpair approx
      CALL dgtapprox(pdl, pd, pdu, er, ei, ncoeff, iseed, berr(i), tol, &
                     i, d, n, td, check)
    ELSE
      !update complex eigenpair approx
      CALL zgtapprox(pdl, pd, pdu, er, ei, ncoeff, iseed, berr(i), tol, &
                     i, d, n, td, check)
    ENDIF
    IF(check) THEN
      IF(it==itmax) PRINT*, i, berr(i)
      EXIT
    ENDIF
  ENDDO
ENDDO
RETURN
END SUBROUTINE dgtlm

!************************************************************************
!			SUBROUTINE DGTAPPROX				*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Laguerre iterate. 				*
!************************************************************************
SUBROUTINE dgtapprox(pdl, pd, pdu, er, ei, ncoeff, iseed, berr, tol, i, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, n, td
REAL(dp), INTENT(IN) :: tol
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: ncoeff(*), pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: j, jmax, jmin
REAL(dp) :: alpha, t
!local arrays
REAL(dp), DIMENSION(n-1) :: adl, adu, c, s, lda
REAL(dp), DIMENSION(n) :: ad
!intrinsic procedures
INTRINSIC :: DABS

!initialize parameters
t=er(i)
!split into 2 cases
IF(DABS(t)>one) THEN
  t=1/t
  CALL drevgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, 0)
  CALL drevseval(ncoeff, DABS(t), alpha, d, 0)
  adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
  DO j=1,n-1
    lda(j)=adl(j)
    IF(DABS(lda(j))<eps) THEN
      lda(j)=eps; adl(j)=eps
    ENDIF
  ENDDO
  CALL dgtqr(adl, ad, adu, c, s, n)
  jmax=DGTJMAX(ad,n)
  jmin=DGTJMIN(ad,n)
  IF(DABS(ad(jmin))<eps) THEN
    !1st stopping criterion met
    berr=DABS(ad(jmin))
    check=.TRUE.
    RETURN
  ENDIF
  CALL dberrapprox(adl, ad, adu, c, s, iseed, berr, n)
  IF(berr<eps) THEN
    !2nd stopping criterion met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criterion met, or it==itmax
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dgtlcorr2(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, &
                 tol, i, d, n, td, check)
  IF(check) THEN
    !3rd stopping criterion met
    RETURN
  ENDIF
ELSE
  CALL dgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, 0)
  CALL dseval(ncoeff, DABS(t), alpha, d, 0)
  adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
  DO j=1,n-1
    lda(j)=adl(j)
    IF(DABS(lda(j))<eps) THEN
      lda(j)=eps; adl(j)=eps
    ENDIF
  ENDDO
  CALL dgtqr(adl, ad, adu, c, s, n)
  jmax=DGTJMAX(ad,n)
  jmin=DGTJMIN(ad,n)
  IF(DABS(ad(jmin))<eps) THEN
    !1st stopping criterion met
    berr=DABS(ad(jmin))
    check=.TRUE.
    RETURN
  ENDIF
  CALL dberrapprox(adl, ad, adu, c, s, iseed, berr, n)
  IF(berr<eps) THEN
    !2nd stopping criterion met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criterion met, or it==itmax
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dgtlcorr1(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, &
                 tol, i, d, n, td, check)
  IF(check) THEN
    !3rd stopping criterion met
    RETURN
  ENDIF
ENDIF

RETURN
END SUBROUTINE dgtapprox

!************************************************************************
!			SUBROUTINE DGTLCORR1				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! general-tridiagonal coeffs of size n and degree d at real number t, 	*
! where |t|<=1. All i-1 previously found eigenvalues are stored in 	*
! (er,ei). The eigenvalue approximation is updated in er(i), ei(i).	*
!************************************************************************
SUBROUTINE dgtlcorr1(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, tol, i, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, n, td
REAL(dp), INTENT(IN) :: alpha, t, tol
!array arguments
REAL(dp), INTENT(IN) :: adl(*), ad(*), adu(*), lda(*), c(*), s(*)
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: k
REAL(dp) :: b0, b1, b2
COMPLEX(dp) :: temp1, temp2, x1, x2, y1, y2
!local arrays
REAL(dp), DIMENSION(n-1) :: bdl, bdu, cdl, cdu
REAL(dp), DIMENSION(n) :: bd, cd, q, v0, v1, v2
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
CALL dgteval(pdl, pd, pdu, t, bdl, bd, bdu, d, n, 1)
CALL dgteval(pdl, pd, pdu, t, cdl, cd, cdu, d, n, 2)
 bdl=bdl/alpha; bd=bd/alpha; bdu=bdu/alpha
 cdl=cdl/alpha; cd=cd/alpha; cdu=cdu/alpha
!compute 1st column of Q^{T}
q(1)=one; q(2:n)=zero
CALL dgtqm('T', c, s, q, n)
!compute v0 and b0 for Hyman's method
b0=ad(n)/q(n)
v0(1:n-1)=b0*q(1:n-1)
v0(n-1)=v0(n-1)-adu(n-1)
IF(n>2) v0(n-2)=v0(n-2)-adl(n-2)
v0(n)=one
CALL dgtrs('N', ad, adu, adl, v0, n-1)
!compute v1 and b1 for Hyman's method
CALL dgtmv('N', bdl, bd, bdu, v0, v1, one, zero, n)
CALL dgtqm('T', c, s, v1, n)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=zero
CALL dgtrs('N', ad, adu, adl, v1, n-1)
!compute b2 for Hyman's method
CALL dgtmv('N', cdl, cd, cdu, v0, v2, one, zero, n)
CALL dgtmv('N', bdl, bd, bdu, v1, v2, two, one, n)
CALL dgtqm('T', c, s, v2, n)
b2=v2(n)/q(n)

!compute temp1 and temp2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=bdl(k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+cdl(k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=p'/p and y2=-(p'/p)'
y1=b1/b0+temp1
y2=y1**2-(b2/b0+2*(b1/b0)*temp1+temp2)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dgtlcorr1

!************************************************************************
!			SUBROUTINE DGTLCORR2				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! general-tridiagonal coeffs of size n and degree d at real number t, 	*
! where |t|>1. All i-1 previously found eigenvalues are stored in 	*
! (er,ei). The eigenvalue approximation is updated in er(i), ei(i).	*
!************************************************************************
SUBROUTINE dgtlcorr2(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, tol, i, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, n, td
REAL(dp), INTENT(IN) :: alpha, t, tol
!array arguments
REAL(dp), INTENT(IN) :: adl(*), ad(*), adu(*), lda(*), c(*), s(*)
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: k
REAL(dp) :: b0, b1, b2
COMPLEX(dp) :: temp1, temp2, x1, x2, y1, y2
!local arrays
REAL(dp), DIMENSION(n-1) :: bdl, bdu, cdl, cdu
REAL(dp), DIMENSION(n) :: bd, cd, q, v0, v1, v2
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
CALL drevgteval(pdl, pd, pdu, t, bdl, bd, bdu, d, n, 1)
CALL drevgteval(pdl, pd, pdu, t, cdl, cd, cdu, d, n, 2)
 bdl=bdl/alpha; bd=bd/alpha; bdu=bdu/alpha
 cdl=cdl/alpha; cd=cd/alpha; cdu=cdu/alpha
!compute 1st column of Q^{T}
q(1)=one; q(2:n)=zero
CALL dgtqm('T', c, s, q, n)
!compute v0 and b0 for Hyman's method
b0=ad(n)/q(n)
v0(1:n-1)=b0*q(1:n-1)
v0(n-1)=v0(n-1)-adu(n-1)
IF(n>2) v0(n-2)=v0(n-2)-adl(n-2)
v0(n)=one
CALL dgtrs('N', ad, adu, adl, v0, n-1)
!compute v1 and b1 for Hyman's method
CALL dgtmv('N', bdl, bd, bdu, v0, v1, one, zero, n)
CALL dgtqm('T', c, s, v1, n)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=zero
CALL dgtrs('N', ad, adu, adl, v1, n-1)
!compute b2 for Hyman's method
CALL dgtmv('N', cdl, cd, cdu, v0, v2, one, zero, n)
CALL dgtmv('N', bdl, bd, bdu, v1, v2, two, one, n)
CALL dgtqm('T', c, s, v2, n)
b2=v2(n)/q(n)

!compute temp1 and temp2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=bdl(k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+cdl(k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=revp'/revp and y2=revp''/revp
y1=b1/b0+temp1
y2=(b2/b0+2*(b1/b0)*temp1+temp2)
!compute y1=p'/p and y2=-(p'/p)'
y2=t**2*(n*d-2*t*y1+t**2*(y1**2-y2))
y1=t*(n*d-t*y1)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dgtlcorr2

!************************************************************************
!			SUBROUTINE DBERRAPPROX				*
!************************************************************************
! Compute an approximation to the backward error associated with the 	*
! smallest eigenvalue of the real matrix a, which is stored in qr	*
! form (dgtqr). Result is stored in berr.				*
!************************************************************************
SUBROUTINE dberrapprox(adl, ad, adu, c, s, iseed, berr, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: adl(*), ad(*), adu(*), c(*), s(*)
!local scalars
INTEGER :: k
REAL(dp) :: bnorm, xnorm
!local arrays
REAL(dp), DIMENSION(3) :: temp
REAL(dp), DIMENSION(n) :: x
!intrinsic procedures
INTRINSIC :: MINVAL
!external procedures
REAL(dp) :: dnrm2
EXTERNAL :: dnrm2

DO k=1,3
  !random real vector
  CALL dlarnv(2, iseed, n, x)
  bnorm=dnrm2(n, x, 1)
  !solve matrix equation
  CALL dgtqm('T', c, s, x, n)
  CALL dgtrs('N', ad, adu, adl, x, n)
  xnorm=dnrm2(n, x, 1)
  temp(k)=bnorm/xnorm
ENDDO
berr=MINVAL(temp)
RETURN
END SUBROUTINE dberrapprox

!****************************************************************
!			SUBROUTINE DREVGTEVAL			*
!****************************************************************
! Evaluate real general-tridiagonal matrix polynomial of size n	*
! and degree d at real number t, where |t|>1. Returns evaluation* 
! in (adl, ad, adu).						*
!****************************************************************
SUBROUTINE drevgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
REAL(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: adl(n-1), ad(n), adu(n-1)
!local scalars
INTEGER :: k

IF(der==0) THEN
  adl=pdl(:,1); ad=pd(:,1); adu=pdu(:,1)
  DO k=2,d+1
    adl=t*adl+pdl(:,k)
    ad=t*ad+pd(:,k)
    adu=t*adu+pdu(:,k)
  ENDDO
ELSEIF(der==1) THEN
  adl=d*pdl(:,1); ad=d*pd(:,1); adu=d*pdu(:,1)
  DO k=2,d
    adl=t*adl+(d-k+1)*pdl(:,k)
    ad=t*ad+(d-k+1)*pd(:,k)
    adu=t*adu+(d-k+1)*pdu(:,k)
  ENDDO
ELSE
  adl=d*(d-1)*pdl(:,1); ad=d*(d-1)*pd(:,1); adu=d*(d-1)*pdu(:,1)
  DO k=2,d-1
    adl=t*adl+(d-k+1)*(d-k)*pdl(:,k)
    ad=t*ad+(d-k+1)*(d-k)*pd(:,k)
    adu=t*adu+(d-k+1)*(d-k)*pdu(:,k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE drevgteval

!****************************************************************
!			SUBROUTINE DGTEVAL			*
!****************************************************************
! Evaluate real general-tridiagonal matrix polynomial of size n	*
! and degree d at complex number t, where |t|<=1. Returns 	*
! evaluation in (adl, ad, adu).					*
!****************************************************************
SUBROUTINE dgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
REAL(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: adl(n-1), ad(n), adu(n-1)
!local scalars
INTEGER :: K

IF(der==0) THEN
  adl=pdl(:,d+1); ad=pd(:,d+1); adu=pdu(:,d+1)
  DO k=d,1,-1
    adl=t*adl+pdl(:,k)
    ad=t*ad+pd(:,k)
    adu=t*adu+pdu(:,k)
  ENDDO
ELSEIF(der==1) THEN
  adl=d*pdl(:,d+1); ad=d*pd(:,d+1); adu=d*pdu(:,d+1)
  DO k=d,2,-1
    adl=t*adl+(k-1)*pdl(:,k)
    ad=t*ad+(k-1)*pd(:,k)
    adu=t*adu+(k-1)*pdu(:,k)
  ENDDO
ELSE
  adl=d*(d-1)*pdl(:,d+1); ad=d*(d-1)*pd(:,d+1); adu=d*(d-1)*pdu(:,d+1)
  DO k=d,3,-1
    adl=t*adl+(k-1)*(k-2)*pdl(:,k)
    ad=t*ad+(k-1)*(k-2)*pd(:,k)
    adu=t*adu+(k-1)*(k-2)*pdu(:,k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE dgteval

!****************************************************************
!			SUBROUTINE DGTSTART			*
!****************************************************************
! Compute the initial estimates of the eigenvalues of a matrix	*
! polynomial of degree d with real general-trdiagonal coeffs of	*
! size n. 							*
!****************************************************************
SUBROUTINE dgtstart(pdl, pd, pdu, er, ei, ncoeff, d, die, dze, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, n
INTEGER, INTENT(INOUT) :: die, dze
!array arguments
REAL(dp), INTENT(IN) :: ncoeff(*), pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: k
!local arrays
LOGICAL, DIMENSION(n) :: pivot
REAL(dp), DIMENSION(n-1) :: adl, adu, c, s
REAL(dp), DIMENSION(n) :: ad
!intrinsic procedures
INTRINSIC :: DABS

!infinite eigenspace
adl=pdl(:,d+1)/ncoeff(d+1)
ad=pd(:,d+1)/ncoeff(d+1)
adu=pdu(:,d+1)/ncoeff(d+1)
CALL dgtqr(adl, ad, adu, c, s, n)
!compute PIVOT, die
CALL diepivot('N', adl, ad, adu, pivot, n)
die=0
DO k=1,n
  IF(.NOT.pivot(k)) die=die+1
ENDDO
!zero eigenspace
adl=pdl(:,1)/ncoeff(1)
ad=pd(:,1)/ncoeff(1)
adu=pdu(:,1)/ncoeff(1)
CALL dgtqr(adl, ad, adu, c, s, n)
!compute PIVOT, dze
CALL diepivot('N', adl, ad, adu, pivot, n)
dze=0
DO k=1,n
  IF(.NOT.pivot(k)) dze=dze+1
ENDDO
!initial estiamtes
CALL dienr(pdl, pd, pdu, er, ei, c, s, d, n)
!store zero and inf eigenvalues
er(1:dze)=zero; ei(1:dze)=zero
er(n*d-die+1:n*d)=big; ei(n*d-die+1:n*d)=zero
RETURN
END SUBROUTINE dgtstart

!************************************************************************
!				SUBROUTINE DIENR        		*
!************************************************************************
! Compute initial estimates for the finite eigenvalues of p using the   *
! numerical range of the matrix polynomial. The matrix a and vectors c	*
! and s are used in dgtqm to provide orthogonal vectors. The roots 	*
! of the scalar polynomial produced are computed using dslm.		*
!************************************************************************
SUBROUTINE dienr(pdl, pd, pdu, er, ei, c, s, d, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, n
!array arguments
REAL(dp), INTENT(IN) :: c(*), s(*), pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: i, j
REAL(dp) :: tr, ti
!local arrays
REAL(dp), DIMENSION(d+1) :: berr, coeff
REAL(dp), DIMENSION(n) :: x, y
!external procedures
REAL(dp) :: ddot
EXTERNAL :: ddot

DO i=1,n
  x=zero; x(i)=one
  CALL dgtqm('N', c, s, x, n)
  DO j=1,d+1
    CALL dgtmv('N', pdl(1,j), pd(1,j), pdu(1,j), x, y, one, zero, n)
    coeff(j)=ddot(n,x,1,y,1)
  ENDDO
  CALL dslm(coeff, er((i-1)*d+1), ei((i-1)*d+1), berr, d)
ENDDO
!sort
DO i=2,n*d
  DO j=i,2,-1
    IF(DCMOD(er(j),ei(j))<DCMOD(er(j-1),ei(j-1))) THEN
      tr=er(j); ti=ei(j)
      er(j)=er(j-1); ei(j)=ei(j-1)
      er(j-1)=tr; ei(j-1)=ti
    ELSE
      EXIT
    ENDIF
  ENDDO
ENDDO
RETURN
END SUBROUTINE dienr

!************************************************************************
!				SUBROUTINE DIEPIVOT			*
!************************************************************************
! Compute pivot columns in qr factorization of tridiagonal matrix 	*
! (DGTQR) stored as (ad, adu, adl). If trans='T', then the pivots of the*
! transpose are found.							*
!************************************************************************
SUBROUTINE diepivot(trans, adl, ad, adu, pivot, n)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=1) :: trans
INTEGER, INTENT(IN) :: n
!array arguments
LOGICAL, INTENT(INOUT) :: pivot(*)
REAL(dp), INTENT(IN) :: adl(*), ad(*), adu(*)
!local scalars
INTEGER :: i, j, k
!local arrays
REAL(dp), DIMENSION(3) :: temp
!intrinsic functions
INTRINSIC :: DABS
!external functions
LOGICAL :: lsame
EXTERNAL :: lsame

!initialize
pivot(1:n)=.FALSE.
!check pivots
k=0
IF(lsame(trans,'N')) THEN
  DO i=1,n
    IF(i<(n-1)) THEN
      temp(1)=ad(i); temp(2)=adu(i); temp(3)=adl(i)
    ELSEIF(i<n) THEN
      temp(1)=ad(i); temp(2)=adu(i)
    ELSE
      temp(1)=ad(i)
    ENDIF
    DO j=MAX(i,k+1),MIN(n,i+2)
      IF(DABS(temp(j-i+1))>=eps) THEN
        pivot(j)=.TRUE.
        k=j
        EXIT
      ENDIF
    ENDDO
  ENDDO
ELSE
  DO j=1,n
    IF(j<(n-1)) THEN
      temp(1)=ad(j); temp(2)=adu(j); temp(3)=adl(j)
    ELSEIF(j<n) THEN
      temp(1)=ad(j); temp(2)=adu(j)
    ELSE
      temp(1)=ad(j)
    ENDIF
    DO i=MAX(j,k+1),MIN(n,j+2)
      IF(DABS(temp(i-j+1))>=eps) THEN
        pivot(j)=.TRUE.
        k=j
      ENDIF
    ENDDO
  ENDDO
ENDIF
RETURN
END SUBROUTINE diepivot

!****************************************************************
!			SUBROUTINE DGTMV			*
!****************************************************************
! Compute the matrix-vector product y=alpha*a*x+beta*y, where a	*
! is a real general-tridiagonal matrix stored in adl, ad, adu, 	*
! x and y are vectors, and alpha and beta are scalars. 		*
!****************************************************************
SUBROUTINE dgtmv(trans, adl, ad, adu, x, y, alpha, beta, n)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=1), INTENT(IN) :: trans
INTEGER, INTENT(IN) :: n
REAL(dp), INTENT(IN) :: alpha, beta
!array arguments
REAL(dp), INTENT(IN) :: adl(*), ad(*), adu(*), x(*)
REAL(dp), INTENT(INOUT) :: y(n)
!local scalars
INTEGER :: j
!local arrays
REAL(dp), DIMENSION(3) :: temp
!external functions
LOGICAL :: lsame
EXTERNAL :: lsame

!initiate y, temp
y=beta*y
temp=zero
!compute product
IF(lsame(trans,'N')) THEN
  temp(1)=ad(1)
  temp(2)=adl(1)
  y(1:2)=y(1:2)+alpha*temp(1:2)*x(1)
  DO j=2,n-1
    temp(1)=adu(j-1)
    temp(2)=ad(j)
    temp(3)=adl(j)
    y(j-1:j+1)=y(j-1:j+1)+alpha*temp*x(j)
  ENDDO
  temp(1)=adu(n-1)
  temp(2)=ad(n)
  y(n-1:n)=y(n-1:n)+alpha*temp(1:2)*x(n)
ELSE
  temp(1)=ad(1)
  temp(2)=adu(1)
  y(1:2)=y(1:2)+alpha*temp(1:2)*x(1)
  DO j=2,n-1
    temp(1)=adl(j-1)
    temp(2)=ad(j)
    temp(3)=adu(j)
    y(j-1:j+1)=y(j-1:j+1)+alpha*temp*x(j)
  ENDDO
  temp(1)=adl(n-1)
  temp(2)=ad(n)
  y(n-1:n)=y(n-1:n)+alpha*temp(1:2)*x(n)
ENDIF
RETURN
END SUBROUTINE dgtmv

!****************************************************************
!			SUBROUTINE DGTRS			*
!****************************************************************
! Solve Rx=b, where R is stored in (rd,rdu,rdu2) by DGTQR.	*
! IF trans='N' then Rx=b, if trans='T' then R^{T}x=b, overwrite	*
! solution in b.						*
!****************************************************************
SUBROUTINE dgtrs(trans, rd, rdu, rdu2, b, n)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=1), INTENT(IN) :: trans
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(IN) :: rd(*), rdu(*), rdu2(*)
REAL(dp), INTENT(INOUT) :: b(*)
!local scalars
INTEGER :: j
!local arrays 
REAL(dp), DIMENSION(2) :: temp
!external functions
LOGICAL :: lsame
EXTERNAL :: lsame

IF(n==1) THEN
  b(1)=b(1)/rd(1)
  RETURN
ENDIF

IF(lsame(trans,'N')) THEN
  DO j=n,3,-1
    !create temp
    temp(1)=rdu2(j-2)
    temp(2)=rdu(j-1)
    !compute b(j), b(j-2:j-1)
    b(j)=b(j)/rd(j)
    b(j-2:j-1)=b(j-2:j-1)-b(j)*temp
  ENDDO
  !j=2
  b(2)=b(2)/rd(2)
  b(1)=b(1)-b(2)*rdu(1)
  !j=1
  b(1)=b(1)/rd(1)
ELSE
  DO j=1,n-2
    !create temp
    temp(1)=rdu(j)
    temp(2)=rdu2(j)
    !compute b(j), b(j+1:j+2)
    b(j)=b(j)/rd(j)
    b(j+1:j+2)=b(j+1:j+2)-b(j)*temp
  ENDDO
  !j=n-1
  b(n-1)=b(n-1)/rd(n-1)
  b(n)=b(n)-b(n-1)*rdu(n-1)
  !j=n
  b(n)=b(n)/rd(n)
ENDIF
RETURN
END SUBROUTINE dgtrs

!****************************************************************
!			SUBROUTINE DGTQM			*
!****************************************************************
! Apply Q as stored in (c,s) by DGTQR to a vector x from the	*
! left. If trans='N' then Qx, if trans='T' then Q^{T}x, result	*
! is stored in x.						*
!****************************************************************
SUBROUTINE dgtqm(trans, c, s, x, n)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=1), INTENT(IN) :: trans
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(IN) :: c(*), s(*)
REAL(dp), INTENT(INOUT) :: x(*)
!local scalars
INTEGER :: j
REAL(dp), DIMENSION(2) :: temp
!external functions
LOGICAL :: lsame
EXTERNAL :: lsame

IF(lsame(trans,'N')) THEN
  DO j=n-1,1,-1
    !compute
    temp(1)=c(j)*x(j)-s(j)*x(j+1)
    temp(2)=s(j)*x(j)+c(j)*x(j+1)
    !update
    x(j:j+1)=temp
  ENDDO
ELSE
  DO j=1,n-1
    !compute
    temp(1)=c(j)*x(j)+s(j)*x(j+1)
    temp(2)=-s(j)*x(j)+c(j)*x(j+1)
    !update
    x(j:j+1)=temp
  ENDDO
ENDIF
RETURN
END SUBROUTINE dgtqm

!****************************************************************
!			SUBROUTINE DGTQR			*
!****************************************************************
! Compute qr factorization of real general-tridiagonal matrix 	*
! stored in (adl, ad, adu). Stores Q in (c,s) and R in  	*
! (adl, ad, adu). Main diagoanl of R is stored in ad, super 	*
! diagonal in adu, and super super diagonal in adl.		*
!****************************************************************
SUBROUTINE dgtqr(adl, ad, adu, c, s, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(INOUT) :: adl(*), ad(*), adu(*), c(*), s(*)
!local scalars
INTEGER :: j
REAL(dp) :: a, b, r

DO j=1,n-1
  !plane rotation
  r=DCMOD(ad(j),adl(j))
  IF(r<eps) THEN
    c(j)=one; s(j)=zero
    ad(j)=zero; adl(j)=zero
  ELSE
    c(j)=ad(j)/r; s(j)=adl(j)/r
    !apply transformation
    ad(j)=r
    a=c(j)*adu(j)+s(j)*ad(j+1)
    b=-s(j)*adu(j)+c(j)*ad(j+1)
    adu(j)=a; ad(j+1)=b
    IF(j<(n-1)) THEN
      adl(j)=s(j)*adu(j+1)
      adu(j+1)=c(j)*adu(j+1)
    ELSE
      adl(j)=zero
    ENDIF
  ENDIF
ENDDO
RETURN
END SUBROUTINE dgtqr

!************************************************************************
!				FUNCTION DGTJMIN			*
!************************************************************************
! Compute j such that |ad(j)| is minimized, or first index such that	* 
! |ad(j)|<eps. ad is real diagonal of a real matrix.			*
!************************************************************************
FUNCTION dgtjmin(ad, n) RESULT(jmin)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(IN) :: ad(*)
!local scalars
INTEGER :: j, jmin
!intrinsic procedures
INTRINSIC :: DABS

jmin=1
DO j=2,n
  IF(DABS(ad(jmin))<eps) RETURN

  IF(DABS(ad(j))<DABS(ad(jmin))) THEN
    jmin=j
  ENDIF
ENDDO
RETURN
END FUNCTION dgtjmin

!************************************************************************
!				FUNCTION DGTJMAX			*
!************************************************************************
! Compute j such that |ad(j)| is minimized, or last index such that	* 
! |ad(j)|<eps. ad is main diagaonal of a real matrix.			*
!************************************************************************
FUNCTION dgtjmax(ad, n) RESULT(jmax)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(IN) :: ad(*)
!local scalars
INTEGER :: j, jmax
!intrinsic procedures
INTRINSIC :: DABS

jmax=n
DO j=n-1,1,-1
  IF(DABS(ad(jmax))<eps) RETURN

  IF(DABS(ad(j))<DABS(ad(jmax))) THEN
    jmax=j
  ENDIF
ENDDO
RETURN
END FUNCTION dgtjmax

!************************************************************************
!			SUBROUTINE ZGTAPPROX				*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Laguerre iterate. 				*
!************************************************************************
SUBROUTINE zgtapprox(pdl, pd, pdu, er, ei, ncoeff, iseed, berr, tol, i, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, n, td
REAL(dp), INTENT(IN) :: tol
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: ncoeff(*), pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: j, jmax, jmin
REAL(dp) :: alpha
COMPLEX(dp) :: t
!local arrays
COMPLEX(dp), DIMENSION(n-1) :: adl, adu, c, s, lda
COMPLEX(dp), DIMENSION(n) :: ad
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS

!initialize parameters
t=DCMPLx(er(i),ei(i))
!split into 2 cases
IF(ZABS(t)>one) THEN
  t=1/t
  CALL zrevgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, 0)
  CALL drevseval(ncoeff, ZABS(t), alpha, d, 0)
  adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
  DO j=1,n-1
    lda(j)=adl(j)
    IF(ZABS(lda(j))<eps) THEN
      lda(j)=eps; adl(j)=eps
    ENDIF
  ENDDO
  CALL zgtqr(adl, ad, adu, c, s, n)
  jmax=ZGTJMAX(ad,n)
  jmin=ZGTJMIN(ad,n)
  IF(ZABS(ad(jmin))<eps) THEN
    !1st stopping criterion met
    berr=ZABS(ad(jmin))
    check=.TRUE.
    RETURN
  ENDIF
  CALL zberrapprox(adl, ad, adu, c, s, iseed, berr, n)
  IF(berr<eps) THEN
    !2nd stopping criterion met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criterion met, or it==itmax
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zgtlcorr2(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, &
                 tol, i, d, n, td, check)
  IF(check) THEN
    !3rd stopping criterion met
    RETURN
  ENDIF
ELSE
  CALL zgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, 0)
  CALL dseval(ncoeff, ZABS(t), alpha, d, 0)
  adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
  DO j=1,n-1
    lda(j)=adl(j)
    IF(ZABS(lda(j))<eps) THEN
      lda(j)=eps; adl(j)=eps
    ENDIF
  ENDDO
  CALL zgtqr(adl, ad, adu, c, s, n)
  jmax=ZGTJMAX(ad,n)
  jmin=ZGTJMIN(ad,n)
  IF(ZABS(ad(jmin))<eps) THEN
    !1st stopping criterion met
    berr=ZABS(ad(jmin))
    check=.TRUE.
    RETURN
  ENDIF
  CALL zberrapprox(adl, ad, adu, c, s, iseed, berr, n)
  IF(berr<eps) THEN
    !2nd stopping criterion met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criterion met, or it==itmax
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zgtlcorr1(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, &
                 tol, i, d, n, td, check)
  IF(check) THEN
    !3rd stopping criterion met
    RETURN
  ENDIF
ENDIF

RETURN
END SUBROUTINE zgtapprox

!************************************************************************
!			SUBROUTINE ZGTLCORR1				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! general-tridiagonal coeffs of size n and degree d at complex number t,*
! where |t|<=1. All i-1 previously found eigenvalues are stored in 	*
! (er,ei). The eigenvalue approximation is updated in er(i), ei(i).	*
!************************************************************************
SUBROUTINE zgtlcorr1(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, tol, i, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, n, td
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(INOUT) :: t
!array arguments
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: adl(*), ad(*), adu(*), lda(*), c(*), s(*)
!local scalars
INTEGER :: k
COMPLEX(dp) :: b0, b1, b2, temp1, temp2, x1, x2, y1, y2
!local arrays
COMPLEX(dp), DIMENSION(n-1) :: bdl, bdu, cdl, cdu
COMPLEX(dp), DIMENSION(n) :: bd, cd, q, v0, v1, v2
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
CALL zgteval(pdl, pd, pdu, t, bdl, bd, bdu, d, n, 1)
CALL zgteval(pdl, pd, pdu, t, cdl, cd, cdu, d, n, 2)
 bdl=bdl/alpha; bd=bd/alpha; bdu=bdu/alpha
 cdl=cdl/alpha; cd=cd/alpha; cdu=cdu/alpha
!compute 1st column of Q^{T}
q(1)=cone; q(2:n)=czero
CALL zgtqm('C', c, s, q, n)
!compute v0 and b0 for Hyman's method
b0=ad(n)/q(n)
v0(1:n-1)=b0*q(1:n-1)
v0(n-1)=v0(n-1)-adu(n-1)
IF(n>2) v0(n-2)=v0(n-2)-adl(n-2)
v0(n)=cone
CALL zgtrs('N', ad, adu, adl, v0, n-1)
!compute v1 and b1 for Hyman's method
CALL zgtmv('N', bdl, bd, bdu, v0, v1, cone, czero, n)
CALL zgtqm('C', c, s, v1, n)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=czero
CALL zgtrs('N', ad, adu, adl, v1, n-1)
!compute b2 for Hyman's method
CALL zgtmv('N', cdl, cd, cdu, v0, v2, cone, czero, n)
CALL zgtmv('N', bdl, bd, bdu, v1, v2, ctwo, cone, n)
CALL zgtqm('C', c, s, v2, n)
b2=v2(n)/q(n)

!compute temp1 and temp2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=bdl(k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+cdl(k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=p'/p and y2=-(p'/p)'
y1=b1/b0+temp1
y2=y1**2-(b2/b0+2*(b1/b0)*temp1+temp2)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=ei(i)-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=ei(i)-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE zgtlcorr1

!************************************************************************
!			SUBROUTINE ZGTLCORR2				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! general-tridiagonal coeffs of size n and degree d at complex number t,*
! where |t|>1. All i-1 previously found eigenvalues are stored in 	*
! (er,ei). The eigenvalue approximation is updated in er(i), ei(i).	*
!************************************************************************
SUBROUTINE zgtlcorr2(pdl, pd, pdu, adl, ad, adu, er, ei, lda, c, s, alpha, t, tol, i, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, n, td
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(INOUT) :: t
!array arguments
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: adl(*), ad(*), adu(*), lda(*), c(*), s(*)
!local scalars
INTEGER :: k
COMPLEX(dp) :: b0, b1, b2, temp1, temp2, x1, x2, y1, y2
!local arrays
COMPLEX(dp), DIMENSION(n-1) :: bdl, bdu, cdl, cdu
COMPLEX(dp), DIMENSION(n) :: bd, cd, q, v0, v1, v2
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
CALL zrevgteval(pdl, pd, pdu, t, bdl, bd, bdu, d, n, 1)
CALL zrevgteval(pdl, pd, pdu, t, cdl, cd, cdu, d, n, 2)
 bdl=bdl/alpha; bd=bd/alpha; bdu=bdu/alpha
 cdl=cdl/alpha; cd=cd/alpha; cdu=cdu/alpha
!compute 1st column of Q^{T}
q(1)=cone; q(2:n)=czero
CALL zgtqm('C', c, s, q, n)
!compute v0 and b0 for Hyman's method
b0=ad(n)/q(n)
v0(1:n-1)=b0*q(1:n-1)
v0(n-1)=v0(n-1)-adu(n-1)
IF(n>2) v0(n-2)=v0(n-2)-adl(n-2)
v0(n)=cone
CALL zgtrs('N', ad, adu, adl, v0, n-1)
!compute v1 and b1 for Hyman's method
CALL zgtmv('N', bdl, bd, bdu, v0, v1, cone, czero, n)
CALL zgtqm('C', c, s, v1, n)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=czero
CALL zgtrs('N', ad, adu, adl, v1, n-1)
!compute b2 for Hyman's method
CALL zgtmv('N', cdl, cd, cdu, v0, v2, cone, czero, n)
CALL zgtmv('N', bdl, bd, bdu, v1, v2, ctwo, cone, n)
CALL zgtqm('C', c, s, v2, n)
b2=v2(n)/q(n)

!compute temp1 and temp2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=bdl(k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+cdl(k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=revp'/revp and y2=revp''/revp
y1=b1/b0+temp1
y2=(b2/b0+2*(b1/b0)*temp1+temp2)
!compute y1=p'/p and y2=-(p'/p)'
y2=t**2*(n*d-2*t*y1+t**2*(y1**2-y2))
y1=t*(n*d-t*y1)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=ei(i)-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=ei(i)-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE zgtlcorr2

!************************************************************************
!			SUBROUTINE ZBERRAPPROX				*
!************************************************************************
! Compute an approximation to the backward error associated with the 	*
! smallest eigenvalue of the complex matrix a, which is stored in qr	*
! form (zgtqr). Result is stored in berr.				*
!************************************************************************
SUBROUTINE zberrapprox(adl, ad, adu, c, s, iseed, berr, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
COMPLEX(dp), INTENT(IN) :: adl(*), ad(*), adu(*), c(*), s(*)
!local scalars
INTEGER :: k
REAL(dp) :: bnorm, xnorm
!local arrays
REAL(dp), DIMENSION(3) :: temp
COMPLEX(dp), DIMENSION(n) :: x
!intrinsic procedures
INTRINSIC :: MINVAL
!external procedures
REAL(dp) :: dznrm2
EXTERNAL :: dznrm2

DO k=1,3
  !random real vector
  CALL zlarnv(2, iseed, n, x)
  bnorm=dznrm2(n, x, 1)
  !solve matrix equation
  CALL zgtqm('C', c, s, x, n)
  CALL zgtrs('N', ad, adu, adl, x, n)
  xnorm=dznrm2(n, x, 1)
  temp(k)=bnorm/xnorm
ENDDO
berr=MINVAL(temp)
RETURN
END SUBROUTINE zberrapprox

!****************************************************************
!			SUBROUTINE ZREVGTEVAL			*
!****************************************************************
! Evaluate real general-tridiagonal matrix polynomial of size n	*
! and degree d at complex number t, where |t|>1. Returns 	*
! evaluation in (adl, ad, adu).					*
!****************************************************************
SUBROUTINE zrevgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
COMPLEX(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
COMPLEX(dp), INTENT(INOUT) :: adl(n-1), ad(n), adu(n-1)
!local scalars
INTEGER :: k

IF(der==0) THEN
  adl=pdl(:,1); ad=pd(:,1); adu=pdu(:,1)
  DO k=2,d+1
    adl=t*adl+pdl(:,k)
    ad=t*ad+pd(:,k)
    adu=t*adu+pdu(:,k)
  ENDDO
ELSEIF(der==1) THEN
  adl=d*pdl(:,1); ad=d*pd(:,1); adu=d*pdu(:,1)
  DO k=2,d
    adl=t*adl+(d-k+1)*pdl(:,k)
    ad=t*ad+(d-k+1)*pd(:,k)
    adu=t*adu+(d-k+1)*pdu(:,k)
  ENDDO
ELSE
  adl=d*(d-1)*pdl(:,1); ad=d*(d-1)*pd(:,1); adu=d*(d-1)*pdu(:,1)
  DO k=2,d-1
    adl=t*adl+(d-k+1)*(d-k)*pdl(:,k)
    ad=t*ad+(d-k+1)*(d-k)*pd(:,k)
    adu=t*adu+(d-k+1)*(d-k)*pdu(:,k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zrevgteval

!****************************************************************
!			SUBROUTINE ZGTEVAL			*
!****************************************************************
! Evaluate real general-tridiagonal matrix polynomial of size n	*
! and degree d at complex number t, where |t|<=1. Returns 	*
! evaluation in (adl, ad, adu).					*
!****************************************************************
SUBROUTINE zgteval(pdl, pd, pdu, t, adl, ad, adu, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
COMPLEX(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: pdl(n-1,*), pd(n,*), pdu(n-1,*)
COMPLEX(dp), INTENT(INOUT) :: adl(n-1), ad(n), adu(n-1)
!local scalars
INTEGER :: K

IF(der==0) THEN
  adl=pdl(:,d+1); ad=pd(:,d+1); adu=pdu(:,d+1)
  DO k=d,1,-1
    adl=t*adl+pdl(:,k)
    ad=t*ad+pd(:,k)
    adu=t*adu+pdu(:,k)
  ENDDO
ELSEIF(der==1) THEN
  adl=d*pdl(:,d+1); ad=d*pd(:,d+1); adu=d*pdu(:,d+1)
  DO k=d,2,-1
    adl=t*adl+(k-1)*pdl(:,k)
    ad=t*ad+(k-1)*pd(:,k)
    adu=t*adu+(k-1)*pdu(:,k)
  ENDDO
ELSE
  adl=d*(d-1)*pdl(:,d+1); ad=d*(d-1)*pd(:,d+1); adu=d*(d-1)*pdu(:,d+1)
  DO k=d,3,-1
    adl=t*adl+(k-1)*(k-2)*pdl(:,k)
    ad=t*ad+(k-1)*(k-2)*pd(:,k)
    adu=t*adu+(k-1)*(k-2)*pdu(:,k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zgteval

!****************************************************************
!			SUBROUTINE ZGTMV			*
!****************************************************************
! Compute the matrix-vector product y=alpha*a*x+beta*y, where a	*
! is a complex general-tridiagonal matrix stored in adl, ad, adu*
! x and y are vectors, and alpha and beta are scalars. If trans=*
! 'N' then a=a, if trans='C' then a=a^{H}.			*
!****************************************************************
SUBROUTINE zgtmv(trans, adl, ad, adu, x, y, alpha, beta, n)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=1), INTENT(IN) :: trans
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(IN) :: alpha, beta
!array arguments
COMPLEX(dp), INTENT(IN) :: adl(*), ad(*), adu(*), x(*)
COMPLEX(dp), INTENT(INOUT) :: y(n)
!local scalars
INTEGER :: j
!local arrays
COMPLEX(dp), DIMENSION(3) :: temp
!intrinsic functions
INTRINSIC :: DCONJG
!external functions
LOGICAL :: lsame
EXTERNAL :: lsame

!initiate y, temp
y=beta*y
temp=czero
!compute product
IF(lsame(trans,'N')) THEN
  temp(1)=ad(1)
  temp(2)=adl(1)
  y(1:2)=y(1:2)+alpha*temp(1:2)*x(1)
  DO j=2,n-1
    temp(1)=adu(j-1)
    temp(2)=ad(j)
    temp(3)=adl(j)
    y(j-1:j+1)=y(j-1:j+1)+alpha*temp*x(j)
  ENDDO
  temp(1)=adu(n-1)
  temp(2)=ad(n)
  y(n-1:n)=y(n-1:n)+alpha*temp(1:2)*x(n)
ELSE
  temp(1)=DCONJG(ad(1))
  temp(2)=DCONJG(adu(1))
  y(1:2)=y(1:2)+alpha*temp(1:2)*x(1)
  DO j=2,n-1
    temp(1)=DCONJG(adl(j-1))
    temp(2)=DCONJG(ad(j))
    temp(3)=DCONJG(adu(j))
    y(j-1:j+1)=y(j-1:j+1)+alpha*temp*x(j)
  ENDDO
  temp(1)=DCONJG(adl(n-1))
  temp(2)=DCONJG(ad(n))
  y(n-1:n)=y(n-1:n)+alpha*temp(1:2)*x(n)
ENDIF
RETURN
END SUBROUTINE zgtmv

!****************************************************************
!			SUBROUTINE ZGTRS			*
!****************************************************************
! Solve Rx=b, where R is stored in (rd,rdu,rdu2) by ZGTQR.	*
! IF trans='N' then Rx=b, if trans='C' then R^{H}x=b, overwrite	*
! solution in b.						*
!****************************************************************
SUBROUTINE zgtrs(trans, rd, rdu, rdu2, b, n)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=1), INTENT(IN) :: trans
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(IN) :: rd(*), rdu(*), rdu2(*)
COMPLEX(dp), INTENT(INOUT) :: b(*)
!local scalars
INTEGER :: j
!local arrays 
COMPLEX(dp), DIMENSION(2) :: temp
!intrinsic functions
INTRINSIC :: DCONJG
!external functions
LOGICAL :: lsame
EXTERNAL :: lsame

IF(lsame(trans,'N')) THEN
  IF(n==1) THEN
    b(1)=b(1)/rd(1)
    RETURN
  ENDIF
  DO j=n,3,-1
    !create temp
    temp(1)=rdu2(j-2)
    temp(2)=rdu(j-1)
    !compute b(j), b(j-2:j-1)
    b(j)=b(j)/rd(j)
    b(j-2:j-1)=b(j-2:j-1)-b(j)*temp
  ENDDO
  !j=2
  b(2)=b(2)/rd(2)
  b(1)=b(1)-b(2)*rdu(1)
  !j=1
  b(1)=b(1)/rd(1)
ELSE
  IF(n==1) THEN
    b(1)=b(1)/DCONJG(rd(1))
    RETURN
  ENDIF
  DO j=1,n-2
    !create temp
    temp(1)=DCONJG(rdu(j))
    temp(2)=DCONJG(rdu2(j))
    !compute b(j), b(j+1:j+2)
    b(j)=b(j)/DCONJG(rd(j))
    b(j+1:j+2)=b(j+1:j+2)-b(j)*temp
  ENDDO
  !j=n-1
  b(n-1)=b(n-1)/DCONJG(rd(n-1))
  b(n)=b(n)-b(n-1)*DCONJG(rdu(n-1))
  !j=n
  b(n)=b(n)/DCONJG(rd(n))
ENDIF
RETURN
END SUBROUTINE zgtrs

!****************************************************************
!			SUBROUTINE ZGTQM			*
!****************************************************************
! Apply Q as stored in (c,s) by ZGTQR to a vector x from the	*
! left. If trans='N' then Qx, if trans='C' then Q^{H}x, result	*
! is stored in x.						*
!****************************************************************
SUBROUTINE zgtqm(trans, c, s, x, n)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=1), INTENT(IN) :: trans
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(IN) :: c(*), s(*)
COMPLEX(dp), INTENT(INOUT) :: x(*)
!local scalars
INTEGER :: j
COMPLEX(dp), DIMENSION(2) :: temp
!intrinsic functions
INTRINSIC :: DCONJG
!external functions
LOGICAL :: lsame
EXTERNAL :: lsame

IF(lsame(trans,'N')) THEN
  DO j=n-1,1,-1
    !compute
    temp(1)=c(j)*x(j)-DCONJG(s(j))*x(j+1)
    temp(2)=s(j)*x(j)+DCONJG(c(j))*x(j+1)
    !update
    x(j:j+1)=temp
  ENDDO
ELSE
  DO j=1,n-1
    !compute
    temp(1)=DCONJG(c(j))*x(j)+DCONJG(s(j))*x(j+1)
    temp(2)=-s(j)*x(j)+c(j)*x(j+1)
    !update
    x(j:j+1)=temp
  ENDDO
ENDIF
RETURN
END SUBROUTINE zgtqm

!****************************************************************
!			SUBROUTINE ZGTQR			*
!****************************************************************
! Compute qr factorization of complex general-tridiagonal matrix*
! stored in (adl, ad, adu). Stores Q in (c,s) and R in  	*
! (adl, ad, adu). Main diagoanl of R is stored in ad, super 	*
! diagonal in adu, and super super diagonal in adl.		*
!****************************************************************
SUBROUTINE zgtqr(adl, ad, adu, c, s, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(INOUT) :: adl(*), ad(*), adu(*), c(*), s(*)
!local scalars
INTEGER :: j
REAL(dp) :: r
COMPLEX(dp) :: a, b
!intrinsic procedures
INTRINSIC :: DCONJG, ZABS

DO j=1,n-1
  !plane rotation
  r=DCMOD(ZABS(ad(j)),ZABS(adl(j)))
  IF(r<eps) THEN
    c(j)=cone; s(j)=czero
    ad(j)=czero; adl(j)=czero
  ELSE
    c(j)=ad(j)/r; s(j)=adl(j)/r
    !apply transformation
    ad(j)=r
    a=DCONJG(c(j))*adu(j)+DCONJG(s(j))*ad(j+1)
    b=-s(j)*adu(j)+c(j)*ad(j+1)
    adu(j)=a; ad(j+1)=b
    IF(j<(n-1)) THEN
      adl(j)=DCONJG(s(j))*adu(j+1)
      adu(j+1)=c(j)*adu(j+1)
    ELSE
      adl(j)=zero
    ENDIF
  ENDIF
ENDDO
RETURN
END SUBROUTINE zgtqr

!************************************************************************
!				FUNCTION ZGTJMIN			*
!************************************************************************
! Compute j such that |ad(j)| is minimized, or first index such that	* 
! |ad(j)|<eps. ad is the main diagonal of a complex valued matrix.	*
!************************************************************************
FUNCTION zgtjmin(ad, n) RESULT(jmin)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(IN) :: ad(*)
!local scalars
INTEGER :: j, jmin
!intrinsic procedures
INTRINSIC :: ZABS

jmin=1
DO j=2,n
  IF(ZABS(ad(jmin))<eps) RETURN

  IF(ZABS(ad(j))<ZABS(ad(jmin))) THEN
    jmin=j
  ENDIF
ENDDO
RETURN
END FUNCTION zgtjmin

!************************************************************************
!				FUNCTION ZGTJMAX			*
!************************************************************************
! Compute j such that |ad(j)| is minimized, or last index such that	* 
! |ad(j)|<eps. ad is the main diagonal of a complex valued matrix. 	*
!************************************************************************
FUNCTION zgtjmax(ad, n) RESULT(jmax)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(IN) :: ad(*)
!local scalars
INTEGER :: j, jmax
!intrinsic procedures
INTRINSIC :: ZABS

jmax=n
DO j=n-1,1,-1
  IF(ZABS(ad(jmax))<eps) RETURN

  IF(ZABS(ad(j))<ZABS(ad(jmax))) THEN
    jmax=j
  ENDIF
ENDDO
RETURN
END FUNCTION zgtjmax

END MODULE dgtlmpep_subroutines2
