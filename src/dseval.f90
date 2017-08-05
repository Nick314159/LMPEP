!************************************************************************
!			        SUBROUTINE DSEVAL				                    *
!************************************************************************
! Evaluate scalar polynomial p with real coeffs of degree d and its	    *
! der=0,1,2 derivatives at real number t, where |t|<=1. Returns		    *
! evaluation in a.							                            *
!************************************************************************
SUBROUTINE dseval(p, t, a, d, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der
REAL(8), INTENT(IN) :: t
REAL(8), INTENT(INOUT) :: a
!array arguments
REAL(8), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(d+1)
  DO k=d,1,-1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(d+1)
  DO k=d,2,-1
    a=t*a+(k-1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(d+1)
  DO k=d,3,-1
    a=t*a+(k-1)*(k-2)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE dseval
