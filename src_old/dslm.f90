!>\author Thomas R. Cameron* and Nikolas I. Steckley**
!>\institution *Davidson College and **Portland State University
!>\date 2017
!>\brief <b> DSLM computes the roots of a real polynomial. </b>
!>\par Purpose:
!>\verbatim
!> DSLM computes the roots of a real polynomial using Laguerre's method. 
!>\endverbatim
!>\param[in] p
!>\verbatim Double precision array of dimension (deg+1), contains polynomial coefficients, ordered from constant to leading. \endverbatim
!>\param[in,out] er
!>\verbatim Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
!>\param[in,out] ei
!>\verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim
!>\param[in,out] berr
!>\verbatim  Double precision array of dimension deg, backward error of each eigenvalue approximation.\endverbatim
!>\param[in] deg
!>\verbatim Integer, degree of the polynomial.\endverbatim
!>\note MEMORY: O(deg), FLOPS: O(deg^2)
!************************************************************************
SUBROUTINE dslm(p, deg, er, ei, berr)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: p(*)
DOUBLE PRECISION, INTENT(INOUT)     :: berr(*), er(*), ei(*)
!local scalars
LOGICAL                             :: check
INTEGER                             :: i, it, nzr, nir, td
DOUBLE PRECISION                    :: tol
!local arrays
DOUBLE PRECISION, DIMENSION(deg+1)  :: alpha
!intrinsic procedures
INTRINSIC                           :: DABS, DBLE, DSQRT, MAXVAL

!infinite roots and true degree
nir=0
alpha=DABS(p(1:deg+1))
tol=MAXVAL(alpha)
DO i=deg+1,2,-1
  IF(alpha(i)<eps*tol) THEN
    er(i-1)=big;ei(i-1)=zero
    nir=nir+1
  ELSE
    EXIT
  ENDIF
ENDDO
td=deg-nir
!deg 0, 1, and 2
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
!deg > 2
  !initial estiamtes
  CALL dsstart(alpha,er,ei,td)
  !zero roots
  nzr=0
  DO i=1,deg
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

