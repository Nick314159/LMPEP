SUBROUTINE ZLASGE( Job, n, m, A, ldA )

! Purpose
! =======

! To normalize the columns of the n by m general complex matrix A.

! Each column is normalized in the 2-norm and/or, optionally, the first
! element of largest absolute value in each column of A is real and positive.

! Arguments
! =========

! Job   - (input) INTEGER
!         On entry: if Job = 0, just normalize each column of A in the
!                   2-norm.
!                   if Job = 1, make the first element of largest absolute
!                   value in each column of A real and positive.
!                   if Job /= 0 and Job /= 1, normalize each column of A
!                   in the 2-norm and make the first element of largest
!                   absolute value real and positive.

! n     - (input) INTEGER
!         On entry: the number of rows of the matrix A.
!         Constraint: n >= 0.

! m     - (input) INTEGER
!         On entry: the number of columns of the matrix A.
!         Constraint: m >= 0.

! A     - (input/output) COMPLEX(KIND=wp) array, dimension (ldA,*)
!         Note: the second dimension of the array A must be at
!               least m.
!         On entry: the n by m matrix A.
!         On exit: the normalized matrix.

! ldA   - (input) INTEGER
!         On entry: the leading dimension of the array A.
!         Constraint: ldA >= n.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER                :: wp = KIND(0.0D0)
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)              :: ldA, n, m
       INTEGER, INTENT (IN)              :: Job
!      .. Array Arguments ..
       COMPLEX (KIND=wp), INTENT (INOUT) :: A(ldA,*)
!      .. Local Scalars ..
       INTEGER                           :: j, k
!      .. External Functions ..
       REAL (KIND=wp), EXTERNAL          :: DZNRM2
       INTEGER, EXTERNAL                 :: IZMAXA
!      .. Intrinsic Functions ..
       INTRINSIC                            REAL
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
IF( (n == 0).OR.(m == 0) )THEN
   RETURN
END IF

IF( JOB == 0 )THEN
   DO j = 1, m
      ! Normalize the columns of A so that norm2(A(1:n,j)) = 1,
      ! j = 1, 2, ..., m
      A(1:n,j) = A(1:n,j)/DZNRM2( n, A(1,j), 1 )
   END DO
ELSE IF( JOB == 1 )THEN
   DO j = 1, m
      ! Make the first element of largest absolute value in each column of A
      ! real and positive,
      k = IZMAXA( n, A(1,j), 1 )
      A(1:n,j) = A(1:n,j)/( A(k,j)/ABS( A(k,j) ) )
      A(k,j) = REAL( A(k,j) )
   END DO
ELSE
   DO j = 1, m
      ! Normalize the columns of A so that norm2(A(1:n,j)) = 1 and the first
      ! element of largest absolute value in the column is real and positive,
      ! j = 1, 2, ..., m
      A(1:n,j) = A(1:n,j)/DZNRM2( n, A(1,j), 1 )
      k = IZMAXA( n, A(1,j), 1 )
      A(1:n,j) = A(1:n,j)/( A(k,j)/ABS( A(k,j) ) )
      A(k,j) = REAL( A(k,j) )
   END DO
END IF

RETURN

END SUBROUTINE ZLASGE
