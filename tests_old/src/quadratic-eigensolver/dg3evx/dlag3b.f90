SUBROUTINE DLAG3B( Job, n, rev, rankA, rankC, kr, indxr, kc, indxc,            &
                   alphar, alphai, beta, p, A, ldA, B, ldB, C, ldC,            &
                   nA, nB, nC, V, ldV, beV, INFO )

! Purpose
! =======

! To find the backward errors of a quadratic eigenvalue problem corresponding
! to the left, or right, eigenvectors given as the columns of the n by 2n
! matrix V.
!
! The quadratic eigenvalue problem is given in the form
!
!   ( (alpha(j)**2)*A + (alpha(j)*beta(j))*B + (beta(j)**2)*C )*x(j) = 0,
!
! or
!
!   (y(j)**H)*( (alpha(j)**2)*A + (alpha(j)*beta(j))*B + (beta(j)**2)*C ) = 0,

! where A, B and C are n by n real matrices, alpha and beta are vectors,
! with beta real and x(j) and y(j) are respectively right and left
! eigenvectors.

! alpha is given by  alpha = alphar + i*alphai, if alphai(j) is positive then
! alphai(j+1) is negative, so that alpha(j) and alpha(j+1) are a complex
! conjugate pair.

! This routine is called by DG3EVX; see that routine for further details.

! Arguments
! =========

! Job    - (input) Character(1)
!          On entry: if Job = 'L', compute the backward errors for the left
!                       eigenvectors,
!                    if Job = 'R', compute the backward errors for the right
!                       eigenvectors.
!          Constraint: Job = 'L', or 'R'.

! n      - (input) INTEGER
!          On entry: the order of the matrices A, B and C.
!          Constraint: n >= 0.

! rev    - (input) LOGICAL
!          On entry: rev is .TRUE. if the reversal problem was solved and
!                    is .FALSE. otherwise.

! rankA  - (input) INTEGER
!          On entry: the rank of the matrix A, as computed by DG3EVX.
!          Constraint: 0 <= rankA <= n.

! rankC  - (input) INTEGER
!          On entry: the rank of the matrix C, as computed by DG3EVX.
!                    When rankA >= rankC then it is assumed that the values
!                    alpha(j), j = rankA+rankC+1, rankA+rankC+2, ..., rankA+n
!                    are zero and
!                    beta(j), j = rankA+n+1, rankA+n+2, ..., 2n
!                    are zero.
!                    When rankA < rankC then it is assumed that the values
!                    alpha(j), j = rankA+rankC+1, rankA+rankC+2, ..., rankC+n
!                    are zero and
!                    beta(j), j = rankC+n+1, rankC+n+2, ..., 2n
!                    are zero.
!          Constraint: 0 <= rankC <= n.

! kr    - (input) INTEGER
!         On entry: the number of columns of V that represent the real vectors,
!                   Vr.
!         Constraint: kr >= 0.

! indxr - (input) INTEGER array, dimension (kr)
!         On entry: the indices of the real eigenvalues.

! kc     - (input) INTEGER
!          On entry: the number of pairs columns of V that represent the complex
!          vectors, Vc.  A column of Vc is given in the form V(j) + i*V(j+1).
!          Constraint: kc >= 0.

! indxc  - (input) INTEGER array, dimension (kc)
!          On entry: the indices of the complex eigenvalues.

! alphar - (input) REAL(KIND=wp) array, dimension (rankA+rankC)
!          On entry: alphar(j), j = 1, 2, ..., rankA+rankC, contains the
!                    real parts of the elements of the vector alpha.  

! alphai - (input) REAL(KIND=wp) array, dimension (rankA+rankC)
!          On entry: alphai(j), j = 1, 2, ..., rankA+rankC, contains the
!                    imaginary parts of the elements of the vector alpha.  

! beta   - (input) REAL(KIND=wp) array, dimension (rankA+rankC)
!          On entry: beta(j), j = 1, 2, ..., rankA+rankC, contains the elements
!                    of the vector beta.

! p      - (input) REAL(KIND=wp) array, dimension (2*n)
!          On entry: p(j), j = 1, 2, ..., 2*n, is a scale factor for the
!                    jth backward error, that divides the 2-norm of the residual
!                    corresponding to the jth eigenvalue/eigenvector to give
!                    the backward error.

! A      - (input) REAL(KIND=wp) array, dimension (ldA,*)
!          Note: the second dimension of the array A must be at
!                least n if nA /= zero.
!          On entry: if nA /= zero, the n by n matrix A.

! ldA    - (input) INTEGER
!          On entry: the leading dimension of the array A.
!          Constraint: ldA >= n if nA /= zero.

! B      - (input) REAL(KIND=wp) array, dimension (ldB,*)
!          Note: the second dimension of the array B must be at
!                least n if nB /= zero.
!          On entry: if nB /= zero, the n by n matrix B.

! ldB    - (input) INTEGER
!          On entry: the leading dimension of the array B.
!          Constraint: ldB >= n if nB /= zero.

! C      - (input) REAL(KIND=wp) array, dimension (ldC,*)
!          Note: the second dimension of the array C must be at
!                least n if nC /= zero.
!          On entry: if nC /= zero, the n by n matrix C.

! ldC    - (input) INTEGER
!          On entry: the leading dimension of the array C.
!          Constraint: ldC >= n if nC /= zero.

! nA     - (input) REAL(KIND=wp)
!          On entry: if nA is zero then the matrix A is taken to be the zero
!                    matrix.

! nB     - (input) REAL(KIND=wp)
!          On entry: if nB is zero then the matrix B is taken to be the zero
!                    matrix.

! nC     - (input) REAL(KIND=wp)
!          On entry: if nC is zero then the matrix C is taken to be the zero
!                    matrix.

! V      - (input/output) REAL(KIND=wp) array, dimension (ldV,*)
!          Note: the second dimension of the array V must be at least 2*n.
!          On entry: the n by 2*n matrix of eigenvectors V.

! ldV    - (input) INTEGER
!          On entry: the leading dimension of the array V.
!          Constraint: ldV >= n.

! beV    - (output) REAL(KIND=wp) array, dimension (2*n)
!          On exit: beV(j), j = 1, 2, ..., 2*n is the backward error for the
!                   jth column of V.

! INFO    - (output) INTEGER 
!           On exit: INFO=0 unless the routine detects an error.
!           INFO > 0
!              Allocation of memory failed.  INFO returns the value of the
!              STAT  flag from the Fortran ALLOCATE statement of the compiler
!              with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER             :: wp = KIND(0.0D0)
       REAL (KIND=wp), PARAMETER      :: ONE = 1.0E+0_wp, ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)           :: kc, kr, ldA, ldB, ldC, ldV, n,        &
                                         rankA, rankC
       INTEGER, INTENT (OUT)          :: INFO
       LOGICAL, INTENT (IN)           :: rev
       CHARACTER (1)                  :: Job
       REAL (KIND=wp), INTENT (IN)    :: nA, nB, nC
!      .. Array Arguments ..
       INTEGER, INTENT (IN)           :: indxc(kc), indxr(kr)
       REAL (KIND=wp), INTENT (IN)    :: alphai(rankA+rankC),                  &
                                         alphar(rankA+rankC),                  &
                                         beta(rankA+rankC), p(2*n)
       REAL (KIND=wp), INTENT (INOUT) :: A(ldA,*), B(ldB, *), C(ldC,*), V(ldV,*)
       REAL (KIND=wp), INTENT (OUT)   :: beV(2*n)
!      .. Local Scalars ..
       INTEGER                        :: k, iqz, iwarn_STAT, j, two_n
       CHARACTER (1)                  :: Trans
!      .. Local Arrays ..
       REAL (KIND=wp), ALLOCATABLE    :: Rr(:,:)
       COMPLEX (KIND=wp), ALLOCATABLE :: Rc(:,:)
!      .. External Functions ..
       REAL (KIND=wp), EXTERNAL       :: DNRM2, DZNRM2
!      .. External Subroutines ..
       EXTERNAL                          DGEMM, DLAG3C, DLAG3R, DLASGE
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( n == 0 )THEN
   RETURN
END IF

two_n = 2*n
iqz = rankA + rankC
IF( (Job == 'R').OR.(job == 'r') )THEN
   Trans = 'N'
ELSE
   Trans = 'T'
END IF

ALLOCATE( Rr(n,kr), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

! The kr backward errors corresponding to the real eigenvalues
CALL DLAG3R( Trans, n, iqz, kr, indxr, alphar, beta, A, ldA, B, ldB, C,ldC,    &
             nA, nB, nC, V, ldV, Rr, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, kr
   k = indxr(j)
   beV(k) = DNRM2( n, Rr(1,j), 1 )/p(k)
END DO

DEALLOCATE( Rr, STAT = iwarn_STAT )

ALLOCATE( Rc(n,kc), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

! The kc backward errors corresponding to the complex eigenvalues
CALL DLAG3C( Trans, n, iqz, kc, indxc, alphar, alphai, beta,                   &
             A, ldA, B, ldB, C,ldC, nA, nB, nC, V, ldV, Rc, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, kc
   k = indxc(j)
   beV(k) = DZNRM2( n, Rc(1,j), 1 )/p(k)
   beV(k+1) = beV(k)
END DO

DEALLOCATE( Rc, STAT = iwarn_STAT )

IF( iqz < two_n )THEN
   ALLOCATE( Rr(n,two_n-iqz), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = iwarn_STAT
      RETURN
   END IF

   ! Backward errors corresponding to the deflated eigenvalues.
   ! For zero eigenvalues, R(:,j) = C*V(:,j) and for infinite eigenvalues,
   ! R(:,j) = A*V(:,j)
   IF( rev )THEN

      CALL DGEMM( Trans, 'NoTrans', n, n-rankA, n, ONE, C, ldC,                &
                  V(1,iqz+1), ldV, ZERO, Rr, n )
      IF( rankC < n )THEN
         CALL DGEMM( Trans, 'NoTrans', n, n-rankC, n, ONE, A, ldA,             &
                     V(1,n+rankC+1), ldV, ZERO, Rr(1,n-rankA+1), n )
      END IF
   ELSE

      IF( rankC < n )THEN
         CALL DGEMM( Trans, 'NoTrans', n, n-rankC, n, ONE, C, ldC,             &
                     V(1,iqz+1), ldV, ZERO, Rr, n )
      END IF
      IF( rankA < n )THEN
         CALL DGEMM( Trans, 'NoTrans', n, n-rankA, n, ONE, A, ldA,             &
                     V(1,n+rankA+1), ldV, ZERO, Rr(1,n-rankC+1), n )
      END IF
   END IF

   DO j = iqz+1, two_n
      beV(j) = DNRM2( n, Rr(1,j-iqz), 1 )/p(j)
   END DO

   DEALLOCATE( Rr, STAT = iwarn_STAT )
END IF

RETURN

END SUBROUTINE DLAG3B
