!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> DZMOD returns moduli of complex number. </b>
!>\par Purpose:
!>\verbatim
!> DZMOD returns moduli of complex number a+bi, while avoiding harmul overflow and underflow. 
!>\endverbatim
!>\param[in] a
!>\verbatim Double precision number, real part. \endverbatim
!>\param[in] b
!>\verbatim Double precision number, imaginary part.\endverbatim
!************************************************************************
DOUBLE PRECISION FUNCTION dzmod(a, b)
IMPLICIT NONE
!scalar arguments
DOUBLE PRECISION, INTENT(IN)    :: a, b
!parameters
DOUBLE PRECISION, PARAMETER     :: zero=0.0D0
DOUBLE PRECISION, PARAMETER     :: eps=EPSILON(zero)
!return scalar
DOUBLE PRECISION                :: r

IF(DABS(a)<eps .AND. DABS(b)<eps) THEN
  r=zero
ELSEIF(DABS(a)<DABS(b)) THEN
  r=DABS(b)*DSQRT(1+(a/b)**2)
ELSE
  r=DABS(a)*DSQRT(1+(b/a)**2)
ENDIF

dzmod=r
RETURN
END FUNCTION dzmod
