!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> ZREVSEVAL evaluates the reversal of a complex scalar polynomial at a complex number. </b>
!>\par Purpose:
!>\verbatim
!> ZREVSEVAL calculates revp(t), revp'(t), or revp''(t), where revp(t)=t^(deg)p(t^(-1)) and t is a nonzero complex number. What derivative is taken is determined by the parameter der and the computation is done via Horner's method. 
!>\endverbatim
!>\param[in] p
!>\verbatim Double complex array of dimension (deg+1), contains polynomial coefficients, ordered from constant to leading. \endverbatim
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
SUBROUTINE zrevseval(p, t, deg, der, a)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)             :: deg, der
DOUBLE COMPLEX, INTENT(IN)      :: t
DOUBLE COMPLEX, INTENT(INOUT)   :: a
!array arguments
DOUBLE COMPLEX, INTENT(IN)      :: p(*)
!local scalars
INTEGER                         :: k

IF(der==0) THEN
!no derivative
  a=p(1)
  DO k=2,deg+1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
!1st derivative 
  a=deg*p(1)
  DO k=2,deg
    a=t*a+(deg-k+1)*p(k)
  ENDDO
ELSE
!2nd derivative
  a=deg*(deg-1)*p(1)
  DO k=2,deg-1
    a=t*a+(deg-k+1)*(deg-k)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zrevseval
