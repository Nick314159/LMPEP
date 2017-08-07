!************************************************************************
!			SUBROUTINE ZSEVAL				*
!************************************************************************
! Evaluate scalar polynomial p with real coeffs of degree d, and its	*
! der=0,1,2 derivatives at complex number t, where |t|<=1. Returns 	*
! evaluation in a.							*
!************************************************************************
SUBROUTINE zseval(p, t, a, d, der)
IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(0.0D0), itmax=50
!scalar arguments
INTEGER, INTENT(IN) :: d, der
COMPLEX(dp), INTENT(IN) :: t
COMPLEX(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
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
END SUBROUTINE zseval
