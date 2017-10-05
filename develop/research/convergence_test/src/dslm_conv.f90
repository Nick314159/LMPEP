!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> DSLM computes the roots of a real polynomial. </b>
!>\par Purpose:
!>\verbatim
!> DSLM computes the roots of a real polynomial using Laguerre's method. 
!>\endverbatim
!>\param[in] p
!>\verbatim Double precision array of dimension (deg+1), contains polynomial coefficients, ordered from constant to leading. \endverbatim
!>\param[in] deg
!>\verbatim Integer, degree of the polynomial.\endverbatim
!>\param[out] er
!>\verbatim Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
!>\param[out] ei
!>\verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim
!>\param[out] berr
!>\verbatim  Double precision array of dimension deg, backward error of each eigenvalue approximation.\endverbatim
!>\note MEMORY: O(deg), FLOPS: O(deg^2)
!************************************************************************
SUBROUTINE dslm_conv(p, deg, er, ei, berr, exacteigs, error)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: p(*)
DOUBLE PRECISION, INTENT(OUT)       :: berr(*), er(*), ei(*) 
DOUBLE COMPLEX, DIMENSION(deg), INTENT(IN)          :: exacteigs(*)
DOUBLE COMPLEX, INTENT(INOUT)       :: error(*)
!local scalars
LOGICAL, DIMENSION(deg)             :: check
INTEGER                             :: i, it, itmax
DOUBLE PRECISION                    :: tol
!local arrays
DOUBLE PRECISION, DIMENSION(deg+1)  :: alpha
!intrinsic procedures
INTRINSIC                           :: dabs, epsilon, dcmplx
!parameters
DOUBLE PRECISION, PARAMETER         :: eps=epsilon(0.0D0)
!external subroutines
EXTERNAL                            :: dsstart, dslcorr, dzslcorr
!external functions
DOUBLE PRECISION                    :: dzmod, dznrm2
EXTERNAL                            :: dzmod, dznrm2

itmax = 60
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
    error(it) = dznrm2(deg, dcmplx(er(1:deg), ei(1:deg))-exacteigs(1:deg), 1)
ENDDO
RETURN
END SUBROUTINE dslm_conv
