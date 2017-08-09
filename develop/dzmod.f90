!************************************************************************
!                           FUNCTION DZMOD      			            *
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!                           Date: 8/9/2017                              *
!************************************************************************
! Compute the moduli of a complex number while avoid harmul overflow    *
! and underflow.                                                        *
!************************************************************************
! Input Variables:                                                      *
!   a: REAL(re8), real part                                             *
!   b: REAL(re8), imaginary part                                        *
!                                                                       *
! Output Variables:                                                     *
!   r: moduli of the complex number a+bi                                *   
!                                                                       *
! Memory: O(1), FLOPS: O(1)                                             *
!************************************************************************
FUNCTION dzmod(a, b) RESULT(r)
USE util
IMPLICIT NONE
!scalar arguments
REAL(KIND=re8), INTENT(IN)  :: a, b
!return scalar
REAL(KIND=re8)              :: r

IF(DABS(a)<eps .AND. DABS(b)<eps) THEN
  r=zero
  RETURN
ENDIF

IF(DABS(a)<DABS(b)) THEN
  r=DABS(b)*DSQRT(1+(a/b)**2)
ELSE
  r=DABS(a)*DSQRT(1+(b/a)**2)
ENDIF
RETURN
END FUNCTION dzmod
