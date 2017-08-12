!>\brief <b>Test the convexity of the angle formed by (IL,A(IL)), (I,A(I)), (IR,A(IR)) at the vertex (I,A(I)), up to within the tolerance TOLER. If convexity holds then the function is set to .TRUE., otherwise CTEST=.FALSE. The parameter TOLER is set to 0.4 by default</b>
!>\note Borrowed from \ref bini
!>\param[in] a
!>\verbatim vector of double \endverbatim
!>\param[in] "IL,I,IR"
!>\verbatim integers such that IL < I < IR \endverbatim
!>\return 
!>\verbatim
!>     .TRUE. if the angle formed by (IL,A(IL)), (I,A(I)), (IR,A(IR)) at 
!>            the vertex (I,A(I)), is convex up to within the tolerance  
!>            TOLER, i.e., if                                            
!>            (A(I)-A(IL))*(IR-I)-(A(IR)-A(I))*(I-IL)>TOLER.             
!>     .FALSE.,  otherwise.
!>\endverbatim
!************************************************************************
FUNCTION ctest(a, il, i, ir) RESULT(OK)
IMPLICIT NONE
INTEGER, INTENT(IN)   :: i, il, ir
DOUBLE PRECISION, INTENT(IN)      :: a(*)
LOGICAL                         :: OK

! Local variables
DOUBLE PRECISION                  :: s1, s2
DOUBLE PRECISION, PARAMETER       :: toler = 0.4D0

s1 = a(i) - a(il)
s2 = a(ir) - a(i)
s1 = s1*(ir-i)
s2 = s2*(i-il)
OK = .FALSE.
IF(s1 > (s2+toler)) OK = .TRUE.
RETURN
END FUNCTION ctest
