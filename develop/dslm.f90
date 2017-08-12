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
!>\param[in,out] er
!>\verbatim Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
!>\param[in,out] ei
!>\verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim
!>\param[in,out] berr
!>\verbatim  Double precision array of dimension deg, backward error of each eigenvalue approximation.\endverbatim
!>\param[in] deg
!>\verbatim Integer, degree of the polynomial.\endverbatim
!>\note MEMORY: O(deg), FLOPS: O(deg^2)
!************************************************************************
SUBROUTINE dslm(p, deg, er, ei, berr)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: p(*)
DOUBLE PRECISION, INTENT(INOUT)     :: berr(*), er(*), ei(*)
!local scalars
LOGICAL                             :: conv
INTEGER                             :: i, it
DOUBLE PRECISION                    :: tol
!local arrays
DOUBLE PRECISION, DIMENSION(deg+1)  :: alpha
!intrinsic procedures
INTRINSIC                           :: dabs, epsilon
!parameters
INTEGER, PARAMETER                  :: itmax=25
DOUBLE PRECISION, PARAMETER         :: eps=epsilon(0.0D0)
!external subroutines
EXTERNAL                            :: dsstart, dslcorr, zslcorr
!external functions
DOUBLE PRECISION                    :: dzmod
EXTERNAL                            :: dzmod

!alpha
DO i=1,deg+1
    alpha(i)=dabs(p(i))
ENDDO
!initial estimates
CALL dsstart(alpha, deg, er, ei)
!Laguerre's method
DO i=1,deg
    conv=.FALSE.
    DO it=1,itmax
        tol=eps*dzmod(er(i),ei(i))
        IF(dabs(ei(i))<tol) THEN
            CALL dslcorr(p, alpha, tol, deg, i, conv, er, ei, berr(i))
        ELSE
            CALL zslcorr(p, alpha, tol, deg, i, conv, er, ei, berr(i))
        ENDIF
        IF(conv) THEN
            EXIT
        ENDIF
    ENDDO
    PRINT*, i, it, berr(i), er(i), ei(i)
ENDDO
RETURN
END SUBROUTINE dslm
