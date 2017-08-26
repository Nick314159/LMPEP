SUBROUTINE DLASGE( Job, n, m, ai, A, ldA, INFO )

! Purpose
! =======

! To normalize the columns of the n by m real matrix A, where column j is
! regarded as real if ai(j) is zero and columns j and j+1 are regarded as
! complex if ai(j) is non-zero, that is y = ( A(:,j) + i*A(:,j+1) ) is taken
! to be a complex vector, and is normalized as such.

! Each column, or double column, is normalized in the 2-norm and/or,
! optionally, the first element of largest absolute value in each column of A
! is positive and real in the complex case.

! Arguments
! =========

! Job   - (input) INTEGER
!         On entry: if Job = 0, just normalize each column of A, or double
!                   column, in the 2-norm.
!                   if Job = 1, make the first element of largest absolute
!                   value in each column of A positive and real in the complex
!                   case.
!                   if Job /= 0 and Job /= 1, normalize each column of A
!                   in the 2-norm and make the first element of largest
!                   absolute value positive and real in the complex case.

! n     - (input) INTEGER
!         On entry: the number of rows of the matrix A.
!         Constraint: n >= 0.

! m     - (input) INTEGER
!         On entry: the number of columns of the matrix A.
!         Constraint: m >= 0.

! ai    - (input) REAL(KIND=wp) array, dimension (m)
!         On entry: ai(j) should be zero if column j of A is to be regarded as
!                   a real vector and should be non-zero if columns j and j+1
!                   of A are to be regarded as a complex vector.  In the complex
!                   case it is immaterial what is contained in ai(j+1).
!                   

! A     - (input/output) REAL(KIND=wp) array, dimension (ldA,*)
!         Note: the second dimension of the array A must be at
!               least m.
!         On entry: the n by m matrix A.
!         On exit: the normalized matrix.

! ldA   - (input) INTEGER
!         On entry: the leading dimension of the array A.
!         Constraint: ldA >= n.

! INFO   - (output) INTEGER 
!          On exit: INFO=0 unless the routine detects an error.
!          INFO > 0
!             Allocation of memory failed.  INFO returns the value of the
!             STAT  flag from the Fortran ALLOCATE statement of the compiler
!             with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER             :: wp = KIND(0.0D0)
       REAL (KIND=wp), PARAMETER      :: ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)           :: ldA, n, m
       INTEGER, INTENT (IN)           :: Job
       INTEGER, INTENT (OUT)          :: INFO
!      .. Array Arguments ..
       REAL (KIND=wp), INTENT (IN)    :: ai(m)
       REAL (KIND=wp), INTENT (INOUT) :: A(ldA,*)
!      .. Local Scalars ..
       INTEGER                        :: iwarn_STAT, j, k
!      .. Local Arrays ..
       COMPLEX (KIND=wp), ALLOCATABLE :: vecc(:)
!      .. External Functions ..
       INTEGER, EXTERNAL              :: IDAMAX, IZMAXA
       REAL (KIND=wp), EXTERNAL       :: DNRM2, DZNRM2
!      .. Intrinsic Functions ..
       INTRINSIC                         CMPLX, REAL
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0
IF( (n == 0).OR.(m == 0) )THEN
   RETURN
END IF

ALLOCATE( vecc(n), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

IF( JOB == 0 )THEN
   ! Normalize the columns of A
   j = 1
   DO
      IF( j > m ) EXIT
      IF( ai(j) == ZERO )THEN
         ! Normalize the real vectors so that norm2(y) = 1
         A(1:n,j) = A(1:n,j)/DNRM2( n, A(1,j), 1 )
         j = j + 1
      ELSE
         ! Normalize the complex vectors so that norm2(y) = 1
         vecc(1:n) = CMPLX( A(1:n,j), A(1:n,j+1), KIND = wp )
         vecc(1:n) = vecc/DZNRM2( n, vecc, 1 )
         A(1:n,j) = REAL( vecc )
         A(1:n,j+1) = AIMAG( vecc )
         j = j + 2
      END IF
   END DO
ELSE IF( JOB == 1 )THEN
   j = 1
   DO
      IF( j > m ) EXIT
      IF( ai(j) == ZERO )THEN
         ! Make the first element of largest absolute value of y positive
         k = IDAMAX( n, A(1,j), 1 )
         IF( A(k,j) < ZERO )THEN
            A(1:n,j) = -A(1:n,j)
         END IF
         j = j + 1
      ELSE
         ! Make the first element of largest absolute value of y real and
         ! positive
         vecc(1:n) = CMPLX( A(1:n,j), A(1:n,j+1), KIND = wp )
         k = IZMAXA( n, vecc, 1 )
         vecc(1:n) = vecc/( vecc(k)/ABS( vecc(k) ) )
         vecc(k) = REAL( vecc(k) )
         A(1:n,j) = REAL( vecc )
         A(1:n,j+1) = AIMAG( vecc )
         j = j + 2
      END IF
   END DO
ELSE
   ! Normalize the columns of A
   j = 1
   DO
      IF( j > m ) EXIT
      IF( ai(j) == ZERO )THEN
         ! Normalize the real vectors so that norm2(y) = 1 and the first
         ! element of largest absolute value is positive
         A(1:n,j) = A(1:n,j)/DNRM2( n, A(1,j), 1 )
         k = IDAMAX( n, A(1,j), 1 )
         IF( A(k,j) < ZERO )THEN
            A(1:n,j) = -A(1:n,j)
         END IF
         j = j + 1
      ELSE
         ! Normalize the complex vectors so that norm2(y) = 1 and the first
         ! element of largest absolute value is real and positive
         vecc(1:n) = CMPLX( A(1:n,j), A(1:n,j+1), KIND = wp )
         vecc(1:n) = vecc/DZNRM2( n, vecc, 1 )
         k = IZMAXA( n, vecc, 1 )
         vecc(1:n) = vecc/( vecc(k)/ABS( vecc(k) ) )
         vecc(k) = REAL( vecc(k) )
         A(1:n,j) = REAL( vecc )
         A(1:n,j+1) = AIMAG( vecc )
         j = j + 2
      END IF
   END DO

END IF
DEALLOCATE( vecc, STAT = iwarn_STAT )

RETURN

END SUBROUTINE DLASGE
