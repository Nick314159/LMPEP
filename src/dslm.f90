IMPLICIT NONE
!parameters
INTEGER, PARAMETER :: dp=KIND(0.0D0), itmax=50
REAL(dp), PARAMETER :: eps=EPSILON(0.0_dp), big=HUGE(0.0_dp)
REAL(dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
COMPLEX(dp), PARAMETER :: czero=DCMPLX(zero), cone=DCMPLX(one)

!************************************************************************
!			SUBROUTINE DSLM					*
!************************************************************************
! Compute the roots of the scalar polynomial p with real coeffs of	*
! degree d using Laguerre's method. The backward error is stored in be	*
! the number of iterations per root is stored in iter, and the roots are* 
! stored in (er,ei).				                        *
!                                                    			*
!************************************************************************
SUBROUTINE dslm(p, er, ei, berr, d)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d
!array arguments
REAL(dp), INTENT(IN) :: p(*)
REAL(dp), INTENT(INOUT) :: berr(*), er(*), ei(*)
!local scalars
LOGICAL :: check
INTEGER :: i, it, nzr, nir, td
REAL(dp) :: tol
!local arrays
REAL(dp), DIMENSION(d+1) :: alpha
!intrinsic procedures
INTRINSIC :: DABS, DBLE, DSQRT, MAXVAL

!infinite roots and true degree
nir=0
alpha=DABS(p(1:d+1))
tol=MAXVAL(alpha)
DO i=d+1,2,-1
  IF(alpha(i)<eps*tol) THEN
    er(i-1)=big;ei(i-1)=zero
    nir=nir+1
  ELSE
    EXIT
  ENDIF
ENDDO
td=d-nir
!degree 1 and 2
IF(td==0) THEN
  PRINT*, 'warning: constant polynomial'
ELSEIF(td==1) THEN
  er(1)=-p(1)/p(2)
  ei(1)=zero
ELSEIF(td==2) THEN
  IF(p(2)**2>=4*p(1)*p(3)) THEN
    er(1)=(-p(2)+DSQRT(p(2)**2-4*p(1)*p(3)))/(2*p(3))
    ei(1)=zero
    er(2)=(-p(2)-DSQRT(p(2)**2-4*p(1)*p(3)))/(2*p(3))
    ei(2)=zero
  ELSE
    er(1)=-p(2)/(2*p(3))
    ei(1)=DSQRT(4*p(1)*p(3)-p(2)**2)/(2*p(3))
    er(2)=er(1)
    ei(2)=-ei(1)
  ENDIF
ELSE
  !initial estiamtes
  CALL dsstart(alpha,er,ei,td)
  !zero roots
  nzr=0
  DO i=1,d
    IF(alpha(i)<eps*tol) THEN
      er(i)=zero;ei(i)=zero
      nzr=nzr+1
    ELSE
      EXIT
    ENDIF
  ENDDO
  !Laguerre's method
  DO i=nzr+1,td
    check=.FALSE.
    DO it=1,itmax
      tol=MAX(eps*DCMOD(er(i),ei(i)), eps)
      IF(DABS(ei(i))<tol) THEN
        CALL dslcorr(p, alpha, er, ei, berr(i), tol, check, td, i)
      ELSE
        CALL zslcorr(p, alpha, er, ei, berr(i), tol, check, td, i)
      ENDIF
      IF(check) THEN
        EXIT
      ENDIF
    ENDDO
  ENDDO
ENDIF
RETURN
END SUBROUTINE dslm

