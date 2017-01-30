SUBROUTINE ZLAQP3( m, n, nA, A, ldA, jpvtA, tauA, rankA, tol, INFO )

! Purpose
! =======

! To find the QR factorization, with column pivoting, of the n by n complex
! matrix A, together with the rank of A.  A is factorized as
!
!    A = Q*R*(P**T),
!
! where Q is an n by n unitary matrix, R is an n by n upper triangular matrix
! of rank(A) and P is a permutation matrix.  The factorization is performed by
! ZGEQP3.

! Arguments
! =========

! m     - (input) INTEGER
!         On entry: the number of rows of the matrix A.
!         Constraint: m >= 0.

! n     - (input) INTEGER
!         On entry: the number of columns of the matrix A.
!         Constraint: n >= 0.

! nA    - (input) REAL(KIND=wp)
!         On entry: usually the norm of the matrix A, but may be the norm of a
!                   matrix containing A, see the argument tol below.  If nA is
!                   zero then the matrix A is taken to be the zero matrix.

! A     - (input/output) COMPLEX(KIND=wp) array, dimension (ldA,*)
!         Note: the second dimension of the array A must be at
!               least n if nA /= zero.
!         On entry: if nA /= zero, the m by n matrix A.
!         On exit: the upper triangle of A contains the MIN(m,n) by n upper
!                  triangular matrix R; the elements below the diagonal,
!                  together with the array TAU, represent the unitary matrix
!                  Q as a product of MIN(m,n) elementary reflectors.

! ldA   - (input) INTEGER
!         On entry: the leading dimension of the array A.
!         Constraint: ldA >= m if nA /= zero.

! jpvtA - (output) INTEGER array, dimension (n)
!         On exit: jpvtA contains the pivot indices representing the
!                  permutation such that the jth column of A*P was the
!                  jpvtA(j)th column of A, j = 1, 2, ..., n.

! tauA  - (output) COMPLEX(kind=wp) array, dimension (MIN(m,n))
!         On exit: the scalar factors of the elementary reflectors.  See
!                  ZGEQP3.

! rankA - (output) INTEGER
!         On exit: the rank of the matrix A.  The rank of A is the value (j-1),
!                  where j is the index of the first diagonal element, r(j,j),
!                  of the upper triangular matrix R, such that
!                  ABS( r(j,j) ) <= tol.

! tol   - (input) REAL(kind=wp)
!         On entry: the tolerance for determining the rank of the matrix A.
!                   If tol = zero, then n*eps*nA is used in place of tol,
!                   where eps is the machine precision as returned by routine
!                   call DLAMCH( 'E' ).  Note that if  tol < zero  then
!                   rankA will be returned as n irrespective of the size of the
!                   elements of R.

! INFO  - (output) INTEGER 
!         On exit: INFO=0 unless the routine detects an error.
!         INFO > 0
!            Allocation of memory failed.  INFO returns the value of the
!            STAT  flag from the Fortran ALLOCATE statement of the compiler
!            with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER                :: wp = KIND(0.0D0)
       COMPLEX (KIND=wp), PARAMETER      :: CZERO = (0.0E+0_wp,0.0E+0_wp)
       REAL (KIND=wp), PARAMETER         :: ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)              :: ldA, m, n
       INTEGER, INTENT (OUT)             :: INFO, rankA
       REAL (KIND=wp), INTENT (IN)       :: nA, tol
!      .. Array Arguments ..
       COMPLEX (KIND=wp), INTENT (INOUT) :: A(ldA,*)
       COMPLEX (KIND=wp), INTENT (OUT)   :: tauA(n)
       INTEGER, INTENT (OUT)             :: jpvtA(n)
!      .. Local Scalars ..
       INTEGER                           :: iwarn_STAT, j, k, lwork
       REAL (KIND=wp)                    :: eps, loctol
!      .. Local Arrays ..
       COMPLEX (KIND=wp), ALLOCATABLE    :: work(:)
       REAL (KIND=wp), ALLOCATABLE       :: rwork(:)
       COMPLEX (KIND=wp)                 :: dummy(1)
       REAL (KIND=wp)                    :: rdummy(1)
!      .. External Functions ..
       REAL (KIND=wp), EXTERNAL          :: DLAMCH
!      .. External Subroutines ..
       EXTERNAL                             ZGEQP3
!      .. Intrinsic Functions ..
       INTRINSIC                            MAX, MIN
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( (m == 0).OR.(n == 0) )THEN
   RETURN
END IF

k = MIN(m,n)
rankA = 0
IF( nA /= ZERO )THEN
   loctol = tol
   IF( loctol == ZERO )THEN
      eps = DLAMCH( 'E' )
      loctol = n*eps*nA
   END IF
   jpvtA(1:n) = 0

   lwork = -1
   CALL ZGEQP3( m, n, A, ldA, jpvtA, tauA, dummy, lwork, rdummy, INFO )
   lwork = dummy(1)
   ALLOCATE( rwork(MAX(1,2*n) ), work(lwork), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = iwarn_STAT
      RETURN
   END IF

   ! [QA,RA,PA] = qr(A);
   CALL ZGEQP3( m, n, A, ldA, jpvtA, tauA, work, lwork, rwork, INFO )

   ! rankA is the rank of R (the upper triangle of A)
   DO j = 1, k
      IF( ABS( a(j,j) ) <= loctol )EXIT
      rankA = rankA + 1
   END DO

   DEALLOCATE( rwork, work, STAT = iwarn_STAT )
ELSE
   jpvtA(1:n) = (/ (j, j = 1, n) /)
   tauA(1:k) = CZERO
   IF( tol < ZERO )THEN
      rankA = n
   END IF
END IF

RETURN

END SUBROUTINE ZLAQP3
