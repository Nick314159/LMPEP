SUBROUTINE ZLAG3V( n, m, alpha, beta, A, ldA, B, ldB, C, ldC, nA, nB, nC,      &
                   V1, ldV1, V2, ldV2, V, ldV, INFO )

! Purpose
! =======

! To find the vectors with smallest backward error of a quadratic eigenvalue
! problem from two sets of n element vectors.  The vectors are the columns of
! the two n by m matrices V1 and V2 and the jth column of each is compared,
! with the column of smallest backward error being stored in the jth column
! of V.
!
! The quadratic eigenvalue problem is given in the form
!
!   ( (alpha(j)**2)*A + (alpha(j)*beta(j))*B + (beta(j)**2)*C )*x(j) = 0,
!
! where A, B and C are n by n complex matrices, alpha and beta are vectors,
! with beta real and x(j) is a vector.

! Arguments
! =========

! n     - (input) INTEGER
!         On entry: the order of the matrices A, B and C.
!         Constraint: n >= 0.

! m     - (input) INTEGER
!         On entry: the number of columns of the matrices V1, V2 and V.
!         Constraint: m >= 0.

! alpha - (input) COMPLEX(KIND=wp) array, dimension (m)
!         On entry: alpha(j), j = 1, 2, ..., m, contains the elements of the
!                   vector alpha.

! beta  - (input) REAL(KIND=wp) array, dimension (m)
!         On entry: beta(j), j = 1, 2, ..., m, contains the elements of the
!                   vector beta.

! A     - (input) COMPLEX(KIND=wp) array, dimension (ldA,*)
!         Note: the second dimension of the array A must be at
!               least n if nA /= zero.
!         On entry: if nA /= zero, the n by n matrix A.

! ldA   - (input) INTEGER
!         On entry: the leading dimension of the array A.
!         Constraint: ldA >= n if nA /= zero.

! B     - (input) COMPLEX(KIND=wp) array, dimension (ldB,*)
!         Note: the second dimension of the array B must be at
!               least n if nB /= zero.
!         On entry: if nB /= zero, the n by n matrix B.

! ldB   - (input) INTEGER
!         On entry: the leading dimension of the array B.
!         Constraint: ldB >= n if nB /= zero.

! C     - (input) COMPLEX(KIND=wp) array, dimension (ldC,*)
!         Note: the second dimension of the array C must be at
!               least n if nC /= zero.
!         On entry: if nC /= zero, the n by n matrix C.

! ldC   - (input) INTEGER
!         On entry: the leading dimension of the array C.
!         Constraint: ldC >= n if nC /= zero.

! nA    - (input) REAL(KIND=wp)
!         On entry: if nA is zero then the matrix A is taken to be the zero
!                   matrix.

! nB    - (input) REAL(KIND=wp)
!         On entry: if nB is zero then the matrix B is taken to be the zero
!                   matrix.

! nC    - (input) REAL(KIND=wp)
!         On entry: if nC is zero then the matrix C is taken to be the zero
!                   matrix.

! V1    - (input/output) COMPLEX(KIND=wp) array, dimension (ldV1,*)
!         Note: the second dimension of the array V must be at least m.
!         On entry: the n by m matrix V.
!         On exit: each column of V1 is normalized in the 2-norm.

! ldV1  - (input) INTEGER
!         On entry: the leading dimension of the array V1.
!         Constraint: ldV1 >= n.

! V2    - (input/output) COMPLEX(KIND=wp) array, dimension (ldV2,*)
!         Note: the second dimension of the array V2 must be at least m.
!         On entry: the n by m matrix V2.
!         On exit: each column of V2 is normalized in the 2-norm.

! ldV2  - (input) INTEGER
!         On entry: the leading dimension of the array V2.
!         Constraint: ldV2 >= n.

! V     - (output) COMPLEX(KIND=wp) array, dimension (ldV,*)
!         Note: the second dimension of the array V must be at least m.
!         On exit: the n by m matrix V.  Each column of V is normalized in the
!                  2-norm and the first element of largest absolute value of
!                  each vector is real and positive

! ldV   - (input) INTEGER
!         On entry: the leading dimension of the array V.
!         Constraint: ldV >= n.

! INFO   - (output) INTEGER 
!          On exit: INFO=0 unless the routine detects an error.
!          INFO > 0
!             Allocation of memory failed.  INFO returns the value of the
!             STAT  flag from the Fortran ALLOCATE statement of the compiler
!             with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER                :: wp = KIND(0.0D0)
       REAL (KIND=wp), PARAMETER         :: ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)              :: ldA, ldB, ldC, ldV, ldV1, ldV2,    &
                                            n, m
       INTEGER, INTENT (OUT)             :: INFO
       REAL (KIND=wp), INTENT (IN)       :: nA, nB, nC
!      .. Array Arguments ..
       COMPLEX (KIND=wp), INTENT (IN)    :: alpha(m)
       COMPLEX (KIND=wp), INTENT (INOUT) :: A(ldA,*), B(ldB, *), C(ldC,*),     &
                                            V1(ldV1,*), V2(ldV2,*)
       COMPLEX (KIND=wp), INTENT (OUT)   :: V(ldV, *)
       REAL (KIND=wp), INTENT (IN)       :: beta(m)
!      .. Local Scalars ..
       INTEGER                           :: iwarn_STAT, j
!      .. Local Arrays ..
       COMPLEX (KIND=wp), ALLOCATABLE    :: R(:,:)
       REAL (KIND=wp), ALLOCATABLE       :: abs1(:), abs2(:), res1(:), res2(:)
!      .. External Subroutines ..
       EXTERNAL                             ZLAG3R, ZLASGE
!      .. Intrinsic Functions ..
       INTRINSIC                            ABS, SUM
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( (n == 0).OR.(m == 0) )THEN
   RETURN
END IF

ALLOCATE( abs1(m), abs2(m), res1(m), res2(m), R(n,m), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

! Need to find the vectors with smallest backward error

! Normalize the vectors so that norm2(y) = 1
CALL ZLASGE( 0, n, m, V1, ldV1 )
CALL ZLASGE( 0, n, m, V2, ldV2 )

DO j = 1, m
   abs1(j) = SUM( ABS( V1(1:n,j) ) )
   abs2(j) = SUM( ABS( V2(1:n,j) ) )
END DO

CALL ZLAG3R( 'NoConj', n, m, alpha, beta, A, n, B, ldB, C, n, nA, nB, nC,      &
             V1, ldV1, R, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, m
   IF( abs1(j) /= ZERO )THEN
      res1(j) = SUM( ABS( R(:,j) ) )/abs1(j)
   ELSE
      res1(j) = ZERO
   END IF
END DO

CALL ZLAG3R( 'NoConj', n, m, alpha, beta, A, n, B, ldB, C, n, nA, nB, nC,      &
             V2, ldV2, R, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, m
   IF( abs2(j) /= ZERO )THEN
      res2(j) = SUM( ABS( R(:,j) ) )/abs2(j)
   ELSE
      res2(j) = ZERO
   END IF
END DO
DO j = 1, m
   IF( ( abs1(j) /= ZERO ).AND.( abs2(j) /= ZERO ) )THEN
      IF( res1(j) <= res2(j) )THEN
         V(1:n,j) = V1(1:n,j)
      ELSE
         V(1:n,j) = V2(1:n,j)
      END IF
   ELSE IF( abs1(j) /= ZERO )THEN
      V(1:n,j) = V1(1:n,j)
   ELSE
      V(1:n,j) = V2(1:n,j)
   END IF
END DO

! Make the first element of largest absolute value of each vector real and
! positive
CALL ZLASGE( 1, n, m, V, ldV )

DEALLOCATE( abs1, abs2, res1, res2, R, STAT = iwarn_STAT )

RETURN

END SUBROUTINE ZLAG3V
