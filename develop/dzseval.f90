!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> DZSEVAL evaluates real scalar polynomial at complex number. </b>
!>\par Purpose:
!>\verbatim
!> DZSEVAL calculates p(t), p'(t), or p''(t), where p is a real scalar polynomial and t is a complex number. What derivative is taken is determined by the parameter der and the computation is done via Horner's method. 
!>\endverbatim
!>\param[in] p
!>\verbatim Double precision array of dimension (deg+1), contains polynomial coefficients, ordered from constant to leading. \endverbatim
!>\param[in] t
!>\verbatim Double complex, number to be evaluated.\endverbatim
!>\param[in] deg
!>\verbatim  Integer, degree of the polynomial.\endverbatim
!>\param[in] der
!>\verbatim  Integer, derivative to be taken (0,1,2).\endverbatim
!>\param[out] a
!>\verbatim  Double complex, return value.\endverbatim
!>\note MEMORY: O(deg), FLOPS: O(deg)
!************************************************************************
SUBROUTINE dzseval(p, t, deg, der, a)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)             :: der, deg
DOUBLE COMPLEX, INTENT(IN)      :: t
DOUBLE COMPLEX, INTENT(OUT)     :: a
!array arguments
DOUBLE PRECISION, INTENT(IN)    :: p(*)
!local scalars
INTEGER                         :: k

IF(der==0) THEN
!no derivative
  a=p(deg+1)
  DO k=deg,1,-1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
!1st derivative
  a=deg*p(deg+1)
  DO k=deg,2,-1
    a=t*a+(k-1)*p(k)
  ENDDO
ELSE
!2nd derivative
  a=deg*(deg-1)*p(deg+1)
  DO k=deg,3,-1
    a=t*a+(k-1)*(k-2)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE dzseval
