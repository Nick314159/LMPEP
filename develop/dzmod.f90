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
!intrinsic functions
INTRINSIC                       :: dabs, dsqrt, epsilon
!parameters
DOUBLE PRECISION, PARAMETER     :: zero=0.0D0
DOUBLE PRECISION, PARAMETER     :: eps=epsilon(zero)

IF(dabs(a)<eps .AND. dabs(b)<eps) THEN
  dzmod=zero
ELSEIF(dabs(a)<dabs(b)) THEN
  dzmod=dabs(b)*dsqrt(1+(a/b)**2)
ELSE
  dzmod=dabs(a)*dsqrt(1+(b/a)**2)
ENDIF
RETURN
END FUNCTION dzmod
