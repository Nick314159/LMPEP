!>\author Thomas R. Cameron* and Nikolas I. Steckley**
!>\institution *Davidson College and **Portland State University
!>\date 2017
!>\brief <b> Returns moduli of complex number. </b>
!>\par Purpose:
!>\verbatim
!> Returns moduli of complex number a+bi, while avoiding harmul overflow and underflow. 
!>\endverbatim
!>\param[in] a
!>\verbatim Double precision number, real part. \endverbatim
!>\param[in] b
!>\verbatim Double precision number, imaginary part.\endverbatim
!************************************************************************
FUNCTION dzmod(a, b) RESULT(r)
USE util
IMPLICIT NONE
!scalar arguments
DOUBLE PRECISION, INTENT(IN)    :: a, b
!return scalar
DOUBLE PRECISION                :: r

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
