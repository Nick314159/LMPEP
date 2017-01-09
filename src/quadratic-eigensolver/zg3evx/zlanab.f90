SUBROUTINE ZLANAB( n, alpha, beta, alphan, betan, absa, INFO )

! Purpose
! =======

! To normalize the eigenvalue pairs (alpha(j), beta(j)), j = 1, 2, ..., n
! to lie on the unit circle.

! Arguments
! =========

! n      - (input) INTEGER
!          On entry: the number of eigenvalues
!          Constraint: n >= 0.

! alpha  - (input) COMPLEX(KIND=wp) array, dimension (n)
!          On entry: alpha(j), j = 1, 2, ..., n, contains the elements of the
!                    vector alpha.

! beta   - (input) COMPLEX(KIND=wp) array, dimension (n)
!          On entry: beta(j), j = 1, 2, ..., n, contains the elements of the
!                    vector beta.  The real part of each element of beta is
!                    assumed to be non-negative and the imaginary part of each
!                    element of beta is assumed to be zero.

! alphan - (output) COMPLEX(KIND=wp) array, dimension (n)
!          On exit: alphan(j), j = 1, 2, ..., n, contains the elements of the
!                    normalized vector alpha.

! betan  - (output) REAL(KIND=wp) array, dimension (n)
!          On exit: betan(j), j = 1, 2, ..., n, contains the elements of the
!                    normalized vector beta and each element is non-negative.

! absa   - (output) REAL(KIND-wp) array, dimension (n)
!          On exit: absa(j), j = 1, 2, ..., n, contains ABS( alphan(j) ).

! INFO   - (output) INTEGER 
!          On exit: INFO=0 unless the routine detects an error.
!          INFO > 0
!             Allocation of memory failed.  INFO returns the value of the
!             STAT  flag from the Fortran ALLOCATE statement of the compiler
!             with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: wp = KIND(0.0D0)
       REAL (KIND=wp), PARAMETER       :: ONE = 1.0E+0_wp, ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)            :: n
       INTEGER, INTENT (OUT)           :: INFO
!      .. Array Arguments ..
       COMPLEX (KIND=wp), INTENT (IN)  :: alpha(n), beta(n)
       COMPLEX (KIND=wp), INTENT (OUT) :: alphan(n)
       REAL (KIND=wp), INTENT (OUT)    :: absa(n), betan(n)
!      .. Local Scalars ..
       INTEGER                         :: iwarn_STAT, j
!      .. Local Arrays ..
       REAL (KIND=wp), ALLOCATABLE     :: absab(:)
!      .. External Functions ..
       REAL (KIND=wp), EXTERNAL        :: DLAPY2
!      .. Intrinsic Functions ..
       INTRINSIC                          ABS
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( n == 0 )THEN
   RETURN
END IF

ALLOCATE( absab(n), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

! Normalize the eigenvalues (alpha, beta) to lie on the unit circle
absa(1:n) = ABS( alpha(1:n) )
DO j = 1, n
   absab(j) = DLAPY2( absa(j), beta(j) )
   IF( absab(j) == ZERO )THEN
      absab(j) = ONE
   END IF
END DO
absa(1:n) = absa(1:n)/absab(1:n)
betan(1:n) = beta(1:n)/absab(1:n)
alphan(1:n) = alpha(1:n)/absab(1:n)

DEALLOCATE( absab, STAT = iwarn_STAT )

RETURN

END SUBROUTINE ZLANAB
