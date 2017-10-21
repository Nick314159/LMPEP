SUBROUTINE dslm_conv(p, deg, er, ei, berr, rr, ri, error)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: p(*)
DOUBLE PRECISION, INTENT(OUT)       :: berr(*),er(*),ei(*),rr(*),ri(*),error(*)
!local scalars
LOGICAL, DIMENSION(deg)             :: check
INTEGER                             :: i, it
DOUBLE PRECISION                    :: tol
!local arrays
DOUBLE PRECISION, DIMENSION(deg+1)  :: alpha
!intrinsic procedures
INTRINSIC                           :: dabs, epsilon, dsqrt, maxval
!parameters
INTEGER, PARAMETER                  :: itmax=60
DOUBLE PRECISION, PARAMETER         :: eps=epsilon(0.0D0)
!external subroutines
EXTERNAL                            :: dsstart, dslcorr, dzslcorr
!external functions
DOUBLE PRECISION                    :: dzmod
EXTERNAL                            :: dzmod

!alpha
DO i=1,deg+1
    alpha(i)=dabs(p(i))
ENDDO
!check
DO i=1,deg
    check(i)=.TRUE.
ENDDO
!initial estimates
CALL dsstart(alpha, deg, er, ei)
! initial error
error(1)=maxval(dsqrt((er(1:deg)-rr(1:deg))**2+(ei(1:deg)-ri(1:deg))**2))
!Laguerre's Method
DO it=1,itmax
    DO i=1,deg
        IF(check(i)) THEN
            tol=eps*dzmod(er(i),ei(i))
            IF(dabs(ei(i))<tol) THEN
                CALL dslcorr(p, alpha, tol, deg, i, check(i), er, ei, berr(i))
            ELSE
                CALL dzslcorr(p, alpha, tol, deg, i, check(i), er, ei, berr(i))
            ENDIF
        ENDIF
    ENDDO
    ! error after it iterations
    error(it+1)=maxval(dsqrt((er(1:deg)-rr(1:deg))**2+(ei(1:deg)-ri(1:deg))**2))
ENDDO
RETURN
END SUBROUTINE dslm_conv
