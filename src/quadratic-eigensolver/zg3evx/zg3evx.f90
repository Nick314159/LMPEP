SUBROUTINE ZG3EVX( scal, jobVL, jobVR, sense, tol, n,                          &
                   A, ldA, B, ldB, C, ldC,                                     &
                   alpha, beta, VL, ldVL, VR, ldVR, s, beVL, beVR,             &
                   iwarn, INFO )

! Purpose
! =======

! To solve the quadratic eigenvalue problem

!    (lambda^2*A + lambda*B + C)*x = 0,

! where A, B and C are n by n complex matrices.

! The routine returns the 2*n eigenvalues, lambda(j), j = 1, 2, ..., 2n and can
! optionally return the corresponding right eigenvectors, x(j) and/or left
! eigenvectors, y(j) as well as the condition numbers of the computed
! eigenvalues and backward errors of the computed right and left
! eigenpairs (corresponding eigenvalues and vectors).  A left eigenvector
! satisfies the equation

!    y^H*(lambda^2*A + lambda*B + C) = 0,

! where y^H is the complex conjugate of y.

! lambda is represented as the pair (alpha, beta), such that
! lambda = alpha/beta.  Note that the computation of alpha/beta may overflow and
! indeed beta may be zero.

! Arguments
! =========

! scal   - (input) INTEGER
!          On entry: determines the form of scaling to be performed on
!                    A, B and C.
!          scal = 0: no scaling
!          scal = 1: (the recommended value.) Fan, Lin and Van Dooren scaling
!                     if norm(B)/sqrt(norm(A)*norm(C)) < 10  and no scaling
!                     otherwise, where norm(Z) is the Frobenius norm of Z.
!          scal = 2: Fan, Lin and Van Dooren scaling.
!          scal = 3: tropical scaling with largest root.
!          scal = 4: tropical scaling with smallest root.
!          Constraint: 0 <= scal <= 4.

! jobVL  - (input) CHARACTER(1)
!          On entry: if jobVL = 'N', do not compute left eigenvectors,
!                    if jobVL = 'V', compute the left eigenvectors.
!                    If sense = 1, 2, 4, 5, 6 or 7, then jobVL must be set to
!                    'V'.
!          Constraint: jobVL = 'N' or 'V'.

! jobVR  - (input) CHARACTER(1)
!          On entry: if jobVR = 'N', do not compute right eigenvectors,
!                    if jobVR = 'V', compute the right eigenvectors.
!                    If sense = 1, 3, 4, 5, 6 or 7, then jobVR must be set to
!                    'V'.
!          Constraint: jobVR = 'N' or 'V'.

! sense  - (input) INTEGER
!          On entry: determines whether, or not, condition numbers and backward
!                    errors are computed.
!          sense = 0: do not compute condition numbers, or backward errors.
!          sense = 1: just compute condition numbers for the eigenvalues.
!          sense = 2: just compute backward errors for the left eigenpairs.
!          sense = 3: just compute backward errors for the right eigenpairs.
!          sense = 4: compute backward errors for the left and right eigenpairs.
!          sense = 5: compute condition numbers for the eigenvalues and
!                     backward errors for the left eigenpairs.
!          sense = 6: compute condition numbers for the eigenvalues and
!                     backward errors for the right eigenpairs.
!          sense = 7: compute condition numbers for the eigenvalues and
!                     backward errors for the left and right eigenpairs.
!          Constraint: 0 <= sense <= 7.

! tol    - (input) REAL(KIND=wp)
!          On entry: tol is used as the tolerance for making decisions on rank
!                    in the deflation procedure.  If tol is zero on entry then
!                    n*eps  is used in place of tol, where eps is the machine
!                    precision as returned by routine call DLAMCH( 'E' ).
!                    A diagonal element of a triangular matrix, R, is regarded
!                    as zero if ABS( r(j,j) ) <= tol*size(X), or n*eps*size(X)
!                    when tol is zero, where size(X) is based on the size of the
!                    absolute values of the elements of the matrix X containing
!                    the matrix R.  See reference [3] for the motivation.
!                    If tol is -1.0 on entry then no deflation is attempted.
!                    The recommended value for tol is zero.

! n      - (input) INTEGER
!          On entry: the order of the matrices A, B and C.
!          Constraint: n >= 0.

! A      - (input/output) COMPLEX(KIND=wp) array, dimension (ldA,*)
!          Note: the second dimension of the array A must be at least n.
!          On entry: the n by n matrix A.
!          On exit: A is used as internal workspace, but if either jobVL = 'V',
!                   or jobVR = 'V', then A is restored on exit.

! ldA    - (input) INTEGER
!          On entry: the leading dimension of the array A.
!          Constraint: ldA >= n.

! B      - (input/output) COMPLEX(KIND=wp) array, dimension (ldB,*)
!          Note: the second dimension of the array B must be at least n.
!          On entry: the n by n matrix B.
!          On exit: B is used as internal workspace, but is restored on exit.

! ldB    - (input) INTEGER
!          On entry: the leading dimension of the array B.
!          Constraint: ldB >= n.

! C      - (input/output) COMPLEX(KIND=wp) array, dimension (ldC,*)
!          Note: the second dimension of the array C must be at least n.
!          On entry: the n by n matrix C.
!          On exit: C is used as internal workspace, but if either jobVL = 'V',
!                   or jobVR = 'V', then C is restored on exit.

! ldC    - (input) INTEGER
!          On entry: the leading dimension of the array C.
!          Constraint: ldC >= n.

! alpha  - (output) COMPLEX(KIND=wp) array, dimension (2*n)
!          On exit: alpha(j), j = 1, 2, ..., 2n contains the first part of the
!                   jth eigenvalue pair (alpha(j), beta(j)) of the quadratic
!                   eigenvalue problem.

! beta   - (output) COMPLEX(KIND=wp) array, dimension (2*n)
!          On exit: beta(j), j = 1, 2, ..., 2n contains the second part of the
!                   jth eigenvalue pair (alpha(j), beta(j)) of the quadratic
!                   eigenvalue problem.  Although beta is declared complex, it
!                   is actually real and non-negative.

! VL     - (output) COMPLEX(KIND=wp) array, dimension (ldVL,*)
!          Note: the second dimension of the array VL must be at least
!                2*n if jobVL = 'V'.
!          On exit: if jobVL = 'V', the left eigenvectors y(j) are stored one
!                   after another in the columns of VL, in the same order as the
!                   corresponding eigenvalues.
!                   Each eigenvector will be normalized with length unity and
!                   with the element of largest modulus real and positive.
!          If jobVL = 'N', VL is not referenced.

! ldVL   - (input) INTEGER
!          On entry: the leading dimension of the array VL.
!          Constraint: ldVL >= n if jobVL = 'V'.

! VR     - (output) COMPLEX(KIND=wp) array, dimension (ldVR,*)
!          Note: the second dimension of the array VR must be at least
!                2*n if jobVR = 'V'.
!          On exit: if jobVR = 'V', the right eigenvectors x(j) are stored one
!                   after another in the columns of VR, in the same order as the
!                   corresponding eigenvalues.
!                   Each eigenvector will be normalized with length unity and
!                   with the element of largest modulus real and positive.
!          If jobVR = 'N', VR is not referenced.

! ldVR   - (input) INTEGER
!          On entry: the leading dimension of the array VR.
!          Constraint: ldVR >= n if jobVR = 'V'.

! s      - (output) REAL(KIND=wp) array, dimension (*)
!          Note: the dimension of the array s must be at least 2*n if
!                sense = 1, 5, 6, or 7.
!          Note also: computing the condition numbers of the eigenvalues
!                     requires that both the left and right eigenvectors be
!                     computed.
!          On exit: if sense = 1, 5, 6, or 7, s(j) contains the condition number
!                   for the jth eigenvalue (large condition numbers imply that
!                   the problem is near one with multiple eigenvalues).
!                   Infinite condition numbers are returned as the Fortran
!                   intrinsic function HUGE(x), where x is of type
!                   REAL(KIND=wp).
!          If sense = 0, 2, 3 or 4, s is not referenced.

! beVL   - (output) REAL(KIND=wp) array, dimension (*)
!          Note: the dimension of the array beVL must be at least 2*n if
!                sense = 2, 4, 5, or 7.
!          On exit: if sense = 2, 4, 5, or 7, beVL(j) contains the backward
!                   error estimate for the computed left eigenpair (lambda(j),
!                   y(j)).
!          If sense = 0, 1, 3, or 6, beVL is not referenced

! beVR   - (output) REAL(KIND=wp) array, dimension (*)
!          Note: the dimension of the array beVR must be at least 2*n if
!                sense = 3, 4, 6, or 7.
!          On exit: if sense = 3, 4, 6, or 7, beVR(j) contains the backward
!                   error estimate for the computed right eigenpair
!                   (lambda(j), x(j)).
!          If sense = 0, 1, 2, or 5, beVR is not referenced

! iwarn  - (output) INTEGER
!          On exit: iwarn will be positive if there are warnings, otherwise
!                   iwarn will be 0.
!          If  INFO = 0  then:
!             If  iwarn = 1  then one, or both, of the matrices A and C is zero.
!                In this case no scaling is performed, even if  scal > 0.
!             If  iwarn = 2  then the matrices A and C are singular, or nearly
!                singular, so the problem is potentially ill-posed.
!             If  iwarn = 3  then both the conditions for  iwarn = 1  and
!                iwarn = 2  above, apply.
!             If  iwarn = 4  then  norm(B) >= 10*sqrt( norm(A)*norm(C) )  and
!                backward stability cannot be guaranteed.
!          If  INFO = 2  then  iwarn  returns the value of INFO from ZGGES.
!          If  INFO = 3  then  iwarn  returns the value of INFO from ZGGEV.
!          If  INFO = -999  then  iwarn  returns the value of the  STAT  flag
!             from the Fortran ALLOCATE statement.

! INFO   - (output) INTEGER 
!          On exit: INFO=0 unless the routine detects an error or a warning has
!                   been flagged.
!
!          INFO = -999
!             Allocation of memory failed.  iwarn returns the value of the
!             STAT  flag from the ALLOCATE statement.
!          INFO < 0 and INFO /= -999
!             If INFO = -i, argument i had an illegal value.
!          INFO = 1
!             The quadratic matrix polynomial is nonregular (singular).
!          INFO = 2
!             The QZ iteration failed in ZGGES.  iwarn returns the value of INFO
!             returned by ZGGES.  This failure is unlikely to happen.
!          INFO = 3
!             The QZ iteration failed in ZGGEV.  iwarn returns the value of INFO
!             returned by ZGGEV.  This failure is unlikely to happen.

! Further Comments
! ================

! The quadratic eigenvalue problem is solved by linearizing the problem and
! solving the resulting 2n by 2n generalized eigenvalue problem.  The chosen
! linearization is combined with eigenvalue parameter scaling to have favourable
! conditioning and backward stability properties.  An initial preprocessing step
! is performed that reveals and deflates the zero and infinite eigenvalues
! contributed by singular leading and trailing matrices.

! Infinite eigenvalues have beta set to zero.  Declaring beta as complex, but
! returning beta as real and non-negative is consistent with the LAPACK
! generalized eigenvalue routines.

! The algorithm is backward stable for problems that are not too heavily damped,
! that is, say, norm(B) <= sqrt( norm(A)*norm(C) ).

! If the problem is heavily damped,then the algorithm is backward stable for
! the larger eigenvalues when  scal = 3  is used, and is backward stable for the
! smaller eigenvalues when  scal = 4  is used.  For complete stability, the
! algorithm should be called twice.

! Further details on the algorithm are given in Reference 3.

! References
! ==========

! [1] H.-Y. Fan, W.-W Lin and P. Van Dooren. Normwise scaling of second order
! polynomial matrices. SIAM J. Matrix Anal. Appl. 26, 1, 252--256, 2004.

! [2] S. Gaubert and M, Sharify. Tropical Scaling of polynomial matrices. Lecture
! Notes in Control and Information Sciences Series, vol. 389, Springer-Verlag,
! Berlin, 291--303, 2009.

! [3] S. Hammarling, C. J. Munro and F. Tisseur. An algorithm for the complete
! solution of quadratic eigenvalue problems. ACM Transaction on Mathematical
! Software, Vol. 39, No. 3, Article 18, April 2013. (MIMS Eprint 2011.86,
! Manchester Institute for Mathematical Sciences, School of Mathematics,
! University of Manchester, Manchester M13 9PL, UK, 2011.
! (http://eprints.ma.man.ac.uk/1690/).)

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER                :: wp = KIND(0.0D0)
       COMPLEX (KIND=wp), PARAMETER      :: CONE = (1.0E+0_wp,0.0E+0_wp)
       COMPLEX (KIND=wp), PARAMETER      :: CZERO = (0.0E+0_wp,0.0E+0_wp)
       REAL (KIND=wp), PARAMETER         :: ONE = 1.0E+0_wp, TWO = 2.0E+0_wp
       REAL (KIND=wp), PARAMETER         :: ZERO = 0.0E+0_wp
       CHARACTER (*), PARAMETER          :: SRNAME = 'ZG3EVX'
!      .. Scalar Arguments ..
       REAL (KIND=wp), INTENT (IN)       :: tol
       INTEGER, INTENT (OUT)             :: INFO, iwarn
       INTEGER, INTENT (IN)              :: ldA, ldB, ldC, ldVL, ldVR, n
       INTEGER, INTENT (IN)              :: scal, sense
       CHARACTER (1), INTENT (IN)        :: jobVR, jobVL
!      .. Array Arguments ..
       COMPLEX (KIND=wp), INTENT (INOUT) :: A(ldA,*), B(ldB, *), C(ldC,*)
       COMPLEX (KIND=wp), INTENT (OUT)   :: alpha(2*n), beta(2*n),             &
                                            VL(ldVL,*), VR(ldVR, *)
       REAL (KIND=wp), INTENT (OUT)      :: beVL(*), beVR(*), s(*)
!      .. Local Scalars ..
       COMPLEX (KIND=wp)                 :: asc
       REAL (KIND=wp)                    :: d, delta, eps, g0, g1, g2, gg,     &
                                            inf, loctol, ltol, nA, nAA, nB,    &
                                            nC, nR3, sa, sb, tau, tmp, toltau
       INTEGER                           :: iqz, iwarn_STAT, j, k, ldVSL,      &
                                            ldVSR, lwork, rankA, rankC,        &
                                            rankR3, sdim, two_n
       LOGICAL                           :: lvec, rev, rvec
!      .. Local Arrays ..
       COMPLEX (KIND=wp), ALLOCATABLE    :: AA(:,:), Acopy(:,:), alphan(:),    &
                                            BB(:,:), Ccopy(:,:), DAX(:,:),     &
                                            PA(:,:), PC(:,:), P3(:,:), Q(:,:), &
                                            QAH(:,:), QC(:,:), Q3H(:,:),       &
                                            R3(:,:), tauA(:), tauC(:),         &
                                            tauR3(:), TEMP(:,:), V(:,:),       &
                                            VSL(:,:), VSR(:,:), work(:)
       COMPLEX (KIND=wp)                 :: dummy(2)
       REAL (KIND=wp), ALLOCATABLE       :: abnu(:,:), absa(:), p(:),          &
                                            rbeta(:), rwork(:)
       REAL (KIND=wp)                    :: rdummy(1)
       INTEGER, ALLOCATABLE              :: ip(:), itemp(:), jpvtA(:),         &
                                            jpvtC(:), jpvtR3(:), jpvtT(:) 
       LOGICAL                           :: ldummy(1)
       LOGICAL, ALLOCATABLE              :: select(:)
!      .. External Functions ..
       REAL (KIND=wp), EXTERNAL          :: DLAMCH, DNRM2, ZLANGE
       LOGICAL, EXTERNAL                 :: ZSELCT
!      .. External Subroutines ..
       EXTERNAL                             XERBLA, ZGEMM, ZGGES, ZGGEV,       &
                                            ZLACPY, ZLAG3V, ZLANAB, ZLAQP3,    &
                                            ZLASGE, ZSWAP, ZTGEVC, ZTRSM,      &
                                            ZTZRZF, ZUNGQR, ZUNMQR, ZUNMRZ
!      .. Intrinsic Functions ..
       INTRINSIC                            ABS, CONJG, DOT_PRODUCT, MAX,      &
                                            REAL, SQRT, TRANSPOSE
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0
iwarn = 0

IF( n == 0 )THEN
   RETURN
END IF

IF( jobVL=='N' .OR. jobVL=='n')THEN
   lvec = .FALSE.
ELSE
   lvec = .TRUE.
END IF
IF( jobVR=='N' .OR. jobVR=='n')THEN
   rvec = .FALSE.
ELSE
   rvec = .TRUE.
END IF

IF( (scal < 0).OR.(scal > 4) )THEN
   INFO = -1
ELSE IF( (jobVL /= 'N').AND.(jobVL /= 'n').AND.                                &
         (jobVL /= 'V').AND.(jobVL /= 'v') )THEN
   INFO = -2
ELSE IF( ( (sense == 1).OR.(sense == 2).OR.(sense == 4).OR.(sense == 5).OR.    &
           (sense == 6).OR.(sense == 7) ).AND.(.NOT.lvec) )THEN
   INFO = -2
ELSE IF( (jobVR /= 'N').AND.(jobVR /= 'n').AND.                                &
         (jobVR /= 'V').AND.(jobVR /= 'v') )THEN
   INFO = -3
ELSE IF( ( (sense == 1).OR.(sense == 3).OR.(sense == 4).OR.(sense == 5).OR.    &
           (sense == 6).OR.(sense == 7) ).AND.(.NOT.rvec) )THEN
   INFO = -3
ELSE IF( (sense < 0).OR.(sense > 7) )THEN
   INFO = -4
ELSE IF( n < 0 )THEN
   INFO = -6
ELSE IF( ldA < n )THEN
   INFO = -8
ELSE IF( ldB < n )THEN
   INFO = -10
ELSE IF( ldC < n )THEN
   INFO = -12
ELSE IF( lvec.AND.(ldVL < n) )THEN
   INFO = -16
ELSE IF( rvec.AND.(ldVR < n) )THEN
   INFO = -18
END IF
IF( INFO < 0 )THEN
   CALL XERBLA( SRNAME, -INFO )
   RETURN
END IF

!%%%%%%%%%%%%%%%%%%%%
! Parameter settings
!%%%%%%%%%%%%%%%%%%%%
inf = DLAMCH( 'O' )
toltau = 10
two_n = 2*n

!%%%%%%%%%%
! Scalings
!%%%%%%%%%%
g2 = ZLANGE( 'F', n, n, A, ldA, rdummy )
g1 = ZLANGE( 'F', n, n, B, ldB, rdummy )
g0 = ZLANGE( 'F', n, n, C, ldC, rdummy )

IF( scal == 0 )THEN
   d = ONE
ELSE
   d = max( g0, g1, g2 )
   IF( d == ZERO )THEN
      ! All three matrices are the zero matrix
      d = ONE
   ELSE
      g2 = g2/d
      g1 = g1/d
      g0 = g0/d
   END IF
END IF
tmp = sqrt(g0*g2)
IF( tmp /= 0 )THEN
   tau = g1/tmp
   IF( tau >= toltau )THEN
      iwarn = 4
   END IF
ELSE
   tau = toltau
END IF
nA = g2
nB = g1
nC = g0

! Eigenvalue parameter scaling
IF( nC == ZERO .OR. nA == ZERO )THEN
   ! No parameter scaling will be performed
   iwarn = 1
   gg = ONE
   delta = ONE
ELSE IF( scal == 0 )THEN
   gg = ONE
   delta = ONE
ELSE IF( tau < toltau .OR. scal == 2 )THEN
   ! Fan Lin and Van Dooren scaling
   gg = sqrt(nC/nA)
   delta = TWO/(nC + nB*gg)
ELSE IF( scal == 3 )THEN
   ! Tropical scaling largest root
   gg = nB/nA
   delta = nA/(nB**2)
ELSE IF( scal == 4 )THEN
   ! Tropical scaling smallest root
   gg = nC/nB
   delta = ONE/nC
ELSE
   ! scal == 1  with tau >= toltau
   gg = ONE
   delta = ONE
END IF

IF( (gg /= ONE).OR.(delta /= ONE).OR.(d /= ONE) )THEN
   sb = gg*delta
   sa = gg*sb

   ! Scale A, B and C
   A(1:n,1:n) = (sa/d)*A(1:n,1:n)
   B(1:n,1:n) = (sb/d)*B(1:n,1:n)
   C(1:n,1:n) = (delta/d)*C(1:n,1:n)

   ! Frobenius norm of scaled A,B,C
   nA = sa*nA
   nB = sb*nB
   nC = delta*nC
ELSE
   sa = ONE
   sb = ONE
ENDIF

loctol = tol
eps = DLAMCH( 'E' )
IF( loctol == ZERO )THEN
   loctol = n*eps
END IF

! If necessary, take a copy of the scaled A and C
IF( lvec.OR.rvec )THEN
   ALLOCATE( Acopy(n,n), Ccopy(n,n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   Acopy(1:n,1:n) = A(1:n,1:n)
   Ccopy(1:n,1:n) = C(1:n,1:n)
END IF

!%%%%%%%%%%%%%%%%%%%
! Rank determination
!%%%%%%%%%%%%%%%%%%%

! Allocate the 2n by 2n matrices AA and BB which will be used in the
! QZ algorithm, along with some other arrays
ALLOCATE( AA(two_n,two_n), BB(two_n,two_n), jpvtA(n), jpvtC(n),                &
          tauA(n), tauC(n), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = -999
   iwarn = iwarn_STAT
   RETURN
END IF

! Store -C in AA(n+1:2*n,1:n) before C gets overwritten (only needed when
! rankA = rankC = n, but we don't know the ranks at this stage).
AA(n+1:two_n,1:n) = -C(1:n,1:n)

! [QA,RA,PA] = qr(A);
ltol = loctol*MAX( ONE, nA )
CALL ZLAQP3( n, n, nA, A, ldA, jpvtA, tauA, rankA, ltol, INFO )
IF( INFO > 0 )THEN
   iwarn = INFO
   INFO = -999
END IF

! [QC,RC,PC] = qr(C);
ltol = loctol*MAX( ONE, nB, nC )
CALL ZLAQP3( n, n, nC, C, ldC, jpvtC, tauC, rankC, ltol, INFO )
IF( INFO > 0 )THEN
   iwarn = INFO
   INFO = -999
END IF

rev = rankC > rankA
IF( rev )THEN
   ! Solve for reversal eigenproblem.
   ALLOCATE( itemp(n), TEMP(n,n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   IF( lvec.OR.rvec )THEN
      TEMP(1:n,1:n) = Acopy
      Acopy(1:n,1:n) = Ccopy
      Ccopy(1:n,1:n) = TEMP
   END IF
   TEMP = A(1:n, 1:n)
   A(1:n, 1:n) = C(1:n, 1:n)
   C(1:n, 1:n) = TEMP
   TEMP(1:n,1) = tauA
   tauA(1:n) = tauC
   tauC(1:n) = TEMP(1:n,1)
   itemp(1:n) = jpvtA
   jpvtA(1:n) = jpvtC
   jpvtC(1:n) = itemp
   itemp(1) = rankA
   rankA = rankC
   rankC = itemp(1)
   tmp = nA
   nA = nC
   nC = tmp

   DEALLOCATE( itemp, TEMP, STAT = iwarn_STAT )
END IF

ALLOCATE( PA(1:n,1:n), PC(1:n,1:n), QAH(n,n), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = -999
   iwarn = iwarn_STAT
   RETURN
END IF

! Form the permutation matrices PA and PC
PA(1:n,1:n) = CZERO
DO j = 1, n
   k = jpvtA(j)
   pa(k,j) = CONE
END DO
PC(1:n,1:n) = CZERO
DO j = 1, n
   k = jpvtC(j)
   pc(k,j) = CONE
END DO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Linearization and deflation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Form QA^H (as QAH)
QAH(1:n,1:n) = CZERO
DO j = 1, n
   QAH(j,j) = CONE
END DO

IF( nA /= ZERO )THEN
   lwork = -1
   CALL ZUNMQR( 'Left', 'ConjTrans', n, n, n, A, ldA, tauA, QAH, n,            &
                dummy, lwork, INFO )
   lwork = dummy(1)
   ALLOCATE( work(lwork), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   CALL ZUNMQR( 'Left', 'ConjTrans', n, n, n, A, ldA, tauA, QAH, n,            &
                work, lwork, INFO )

   DEALLOCATE( work, STAT = iwarn_STAT )

   ! Form -QA^H*B in the top n by n part of AA
   IF( nB /= ZERO )THEN

      CALL ZGEMM( 'NoTrans', 'NoTrans', n, n, n, -CONE, QAH, n, B, ldB,        &
                  CZERO, AA, two_n )

   ELSE
      AA(1:n,1:n) = CZERO
   END IF
ELSE
   AA(1:n,1:n) = -B(1:n,1:n)
END IF

IF( rankC == n .AND. rankA == n )THEN
   ! No deflation

   ! Form AA = ( -QA^H*B*PA  QA^H )  and  BB = ( RA  0 )
   !           (    -C*PA     0   )            (  0  I )

   ! Form ( -QA^H*B )*PA
   !      (    -C   )

   ! AA(n+1:2*n,1:n) already contains -C

   !  AA(1:2*n,1:n) = MATMUL( AA(1:2*n,1:n), PA )
   AA(1:two_n,1:n) = AA(1:two_n,jpvtA(:))

   AA(1:n,n+1:two_n) = QAH
   AA(n+1:two_n,n+1:two_n) = CZERO

   BB(1:two_n,1:two_n) = CZERO
   CALL ZLACPY( 'Upper', n, n, A, ldA, BB, two_n )
   DO j = n+1, two_n
      bb(j,j) = ONE
   END DO

   ALLOCATE( V(two_n,two_n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! Form  the orthonormal matrix V = ( PA  0 )
   !                                  (  0  I )
   V(1:n,1:n) = PA
   V(n+1:two_n,1:n) = CZERO
   V(1:two_n,n+1:two_n) = CZERO
   DO j = n+1, two_n
      v(j,j) = CONE
   END DO

   IF( lvec )THEN
      ALLOCATE( Q(two_n,two_n), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! Form the orthonormal matrix  Q = ( QA^H  0 )
      !                                  (  0    I )
      Q(1:n,1:n) = QAH
      Q(n+1:two_n,1:n) = CZERO
      Q(1:two_n,n+1:two_n) = CZERO
      DO j = n+1, two_n
         q(j,j) = CONE
      END DO
   END IF
ELSE IF( rankC < rankA .AND. rankA == n )THEN
   ! Form AA = (  -QA^H*B*PA  QA^H*QC )  and  BB = ( RA  0 )
   !           ( -RC*PC^T*PA        0 )            (  0  I )
   !           (           0        0 )
   !
   ! AA(1:n,1:n) already contains -QA^H*B.  Copy RC to AA and form -RC*PC^T
   AA(n+1:two_n,1:two_n) = CZERO
   CALL ZLACPY( 'Upper', rankC, n, C, ldC, AA(n+1,1), two_n )
   AA(n+1:n+rankC,1:n) = -AA(n+1:n+rankC,1:n)

   ALLOCATE( jpvtT(n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   !  AA(n+1:n+rankC,1:n) = MATMUL( AA(n+1:n+rankC,1:n), TRANSPOSE( PC ) )
   jpvtT( jpvtC(:) ) = (/ (j, j = 1, n) /)
   AA(n+1:n+rankC,1:n) = AA(n+1:n+rankC,jpvtT(:))

   DEALLOCATE( jpvtT, STAT = iwarn_STAT )

   ! Form (  -QA^H*B )*PA
   !      ( -RC*PC^T )

   !  AA(1:n+rankC,1:n) = MATMUL( AA(1:n+rankC,1:n), PA )
   AA(1:n+rankC,1:n) = AA(1:n+rankC,jpvtA(:))

   ! Form QA^H*QC
   AA(1:n,n+1:two_n) = QAH
   IF( nC /= ZERO )THEN
      lwork = -1
      CALL ZUNMQR( 'Right', 'NoConjTrans', n, n, n, C, ldC, tauC,              &
                   AA(1,n+1), two_n, dummy, lwork, INFO )
      lwork = dummy(1)
      ALLOCATE( work(lwork), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      CALL ZUNMQR( 'Right', 'NoConjTrans', n, n, n, C, ldC, tauC,              &
                AA(1,n+1), two_n, work, lwork, INFO )

      DEALLOCATE( work, STAT = iwarn_STAT )
   END IF

   BB(1:two_n,1:two_n) = CZERO
   CALL ZLACPY( 'Upper', n, n, A, ldA, BB, two_n )
   DO j = n+1, two_n
      bb(j,j) = ONE
   END DO

   ALLOCATE( V(two_n,two_n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! Form the orthonormal matrix  V = ( PA  0  )
   !                                  (  0  QC )
   V(1:n,1:n) = PA
   V(n+1:two_n,1:n) = CZERO
   IF( nC /= ZERO )THEN
      V(1:n,n+1:two_n) = CZERO
      V(n+1:two_n,n+1:two_n) = C(1:n,1:n)

      lwork = -1
      CALL ZUNGQR( n, n, n, V(n+1,n+1), two_n, tauC, dummy, lwork, INFO )
      lwork = dummy(1)
      ALLOCATE( work(lwork), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! Form QC in V(n+1:2*n,n+1:2*n)
      CALL ZUNGQR( n, n, n, V(n+1,n+1), two_n, tauC, work, lwork, INFO )

      DEALLOCATE( work, STAT = iwarn_STAT )
   ELSE
      V(1:two_n,n+1:two_n) = CZERO
      DO j = n+1, two_n
         V(j,j) = CONE
      END DO
   END IF

   IF( lvec )THEN
      ALLOCATE( Q(two_n,two_n), QC(n,n-rankC), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! Form the orthonormal matrix  Q = ( QA^H  0   ) and retrieve QC
      !                                  (  0   QC^H )
      Q(1:n,1:n) = QAH
      Q(n+1:two_n,1:n) = CZERO
      Q(1:n,n+1:two_n) = CZERO
      QC(1:n,1:n-rankC) = V(n+1:two_n,n+rankC+1:two_n)
      IF( nC /= ZERO )THEN
         Q(n+1:two_n,n+1:two_n) = TRANSPOSE( CONJG( V(n+1:two_n,n+1:two_n) ) )
      ELSE
         Q(n+1:two_n,n+1:two_n) = V(n+1:two_n,n+1:two_n)
      END IF
   END IF
ELSE
   ! rankC <= rankA < n
   ! A and C are singular, or nearly singular, so the problem is potentially
   ! ill-posed
   iwarn = iwarn + 2

   ! Form the A and B matrices in Qtilde*C2(lambda)*Vtilde (in AA and BB)

   ! Form AA = ( -QA^H*B   QA^H*QC )
   !           ( -RC*PC^T      0   )
   !           (    0          0   )

   ! AA(1:n,1:n) already contains -QA^H*B.  Copy RC to AA and form -RC*PC^T
   AA(n+1:n+rankC,1:n) = CZERO

   CALL ZLACPY( 'Upper', rankC, n, C, ldC, AA(n+1,1), two_n )

   ALLOCATE( jpvtT(n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   !  AA(n+1:n+rankC,1:n) = -AA(n+1:n+rankC,1:n)
   !  AA(n+1:n+rankC,1:n) = MATMUL( AA(n+1:n+rankC,1:n), TRANSPOSE( PC ) )
   jpvtT( jpvtC(:) ) = (/ (j, j = 1, n) /)
   AA(n+1:n+rankC,1:n) = -AA(n+1:n+rankC,jpvtT(:))

   ! Form QA^H*QC and if left eigenvectors are required also form
   ! QC(:,rankC+1:n)
   AA(1:n,n+1:two_n) = QAH

   IF( lvec )THEN
      ALLOCATE( QC(n,n-rankC), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      QC(1:n,1:n-rankC) = CZERO
      DO j = 1, n-rankC
         QC(j+rankC,j) = CONE
      END DO

      IF( nC /= ZERO )THEN
         lwork = -1
         CALL ZUNMQR( 'Right', 'NoConjTrans', n, n, n, C, ldC, tauC,           &
                      AA(1,n+1), two_n, dummy, lwork, INFO )
         CALL ZUNMQR( 'Left', 'NoConjTrans', n, n-rankC, n, C, ldC, tauC,      &
                      QC, n, dummy(2), lwork, INFO )
         lwork = MAX( REAL( dummy(1) ), REAL( dummy(2) ) )
         ALLOCATE( work(lwork), STAT = iwarn_STAT )
         IF( iwarn_STAT /= 0 )THEN
            INFO = -999
            iwarn = iwarn_STAT
            RETURN
         END IF

         CALL ZUNMQR( 'Right', 'NoConjTrans', n, n, n, C, ldC, tauC,           &
                      AA(1,n+1), two_n, work, lwork, INFO )
         CALL ZUNMQR( 'Left', 'NoConjTrans', n, n-rankC, n, C, ldC, tauC,      &
                      QC, n, work, lwork, INFO )

         DEALLOCATE( work, STAT = iwarn_STAT )
      END IF
   ELSE IF( nC /= ZERO )THEN
      lwork = -1
      CALL ZUNMQR( 'Right', 'NoConjTrans', n, n, n, C, ldC, tauC,              &
                   AA(1,n+1), two_n, dummy, lwork, INFO )
      lwork = dummy(1)
      ALLOCATE( work(lwork), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      CALL ZUNMQR( 'Right', 'NoConjTrans', n, n, n, C, ldC, tauC,              &
                   AA(1,n+1), two_n, work, lwork, INFO )

      DEALLOCATE( work, STAT = iwarn_STAT )
   END IF

   AA(n+1:n+rankC,n+1:two_n) = CZERO
   AA(n+rankC+1:two_n,1:two_n) = CZERO

   ! Form BB = ( RA*PA^T  0 )
   !           (   0      0 )
   !           (   0      I )
   BB(1:two_n,1:two_n) = CZERO

   CALL ZLACPY( 'Upper', rankA, n, A, ldA, BB, two_n )

   ! BB(1:rankA,1:n) = MATMUL( BB(1:rankA,1:n), TRANSPOSE( PA ) )
   jpvtT( jpvtA(:) ) = (/ (j, j = 1, n) /)
   BB(1:rankA,1:n) = BB(1:rankA,jpvtT(:))
   DO j = n+1, two_n
      bb(j,j) = ONE
   END DO

   DEALLOCATE( jpvtT, STAT = iwarn_STAT )

   ALLOCATE( jpvtR3(n+rankC), R3(n-rankA,n+rankC), tauR3(n+rankC),             &
             STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! Copy AA(rankA+1:n,1:n+rankC) to R3 and perform a QR factorization of R3
   CALL ZLACPY( 'General', n-rankA, n+rankC, AA(rankA+1,1), two_n, R3, n-rankA )

   ! Compute the Frobenius norm of R3
   nR3 = ZLANGE( 'F', n-rankA, n+rankC, R3, n-rankA, rdummy )
   ! Compute the Frobenius norm of AA
   nAA = ZLANGE( 'F', two_n, two_n, AA, two_n, rdummy )

   ! [Q3,R3,P3] = qr(AA(r2+1:n,1:n+r0));
   ltol = loctol*nAA
   CALL ZLAQP3( n-rankA, n+rankC, nR3, R3, n-rankA, jpvtR3, tauR3, rankR3,     &
                ltol, INFO )
   IF( INFO > 0 )THEN
      iwarn = INFO
      INFO = -999
   END IF

   ALLOCATE( P3(n+rankC,n+rankC), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! Form the permutation matrix P3
   P3(1:n+rankC,1:n+rankC) = CZERO
   DO j = 1, n+rankC
      k = jpvtR3(j)
      p3(k,j) = CONE
   END DO

   DEALLOCATE( jpvtR3, STAT = iwarn_STAT )

   ALLOCATE( Q3H(n-rankA,n-rankA), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! Form Q3^H (as Q3H)
   Q3H(1:n-rankA,1:n-rankA) = CZERO
   DO j = 1, n-rankA
      Q3H(j,j) = CONE
   END DO

   lwork = -1
   CALL ZUNMQR( 'Left', 'ConjTrans', n-rankA, n-rankA, n-rankA, R3, n-rankA,   &
                tauR3, Q3H, n-rankA, dummy, lwork, INFO )
   lwork = dummy(1)
   ALLOCATE( work(lwork), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   CALL ZUNMQR( 'Left', 'ConjTrans', n-rankA, n-rankA, n-rankA, R3, n-rankA,   &
                tauR3, Q3H, n-rankA, work, lwork, INFO )

   DEALLOCATE( work, STAT = iwarn_STAT )

   ! Test for non-regularity
   IF( rankR3 < n-rankA )THEN
      INFO = 1
      RETURN
   END IF

   ! If rankR3 /= n+rankC, reduce R3 to upper triangular form. (This gives the
   ! complete orthogonal decomposition
   IF( rankR3 /= n+rankC )THEN

      lwork = -1
      CALL ZTZRZF( rankR3, n+rankC, R3, n-rankA, tauR3, dummy, lwork, INFO )
      CALL ZUNMRZ( 'Right', 'ConjTrans', n+rankC, n+rankC, rankR3,             &
                   rankC+rankA, R3, n-rankA, tauR3, P3, n+rankC,               &
                   dummy(2), lwork, INFO ) 
      lwork = MAX( REAL( dummy(1) ), REAL( dummy(2) ) )
      ALLOCATE( work(lwork), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      CALL ZTZRZF( rankR3, n+rankC, R3, n-rankA, tauR3, work, lwork, INFO )

      CALL ZUNMRZ( 'Right', 'ConjTrans', n+rankC, n+rankC, rankR3,             &
                   rankC+rankA, R3, n-rankA, tauR3, P3, n+rankC, work, lwork,  &
                   INFO )

      DEALLOCATE( work, STAT = iwarn_STAT )
   ELSE
      P3(1:n+rankC,1:n+rankC) = CZERO
      DO j = 1, n+rankC
        P3(j,j) = CONE
      END DO
   END IF

   ! Now form the matrices A and B of Q*C2(lambda)*V (in AA and BB)

   ! AA(rankA+1:n,1:n+rankC) = [R3 zeros(n-rankA,rankC+rankA)];

   AA(rankA+1:n,1:n+rankC) = CZERO

   CALL ZLACPY( 'Upper', rankR3, rankR3, R3, n-rankA, AA(rankA+1, 1), two_n )

   DEALLOCATE( R3, tauR3, STAT = iwarn_STAT )

   ALLOCATE( TEMP(two_n,n+rankC), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! AA(rankA+1:n,n+rankC+1:2*n) = Q3^H*AA(rankA+1:n,n+rankC+1:2*n);
   ! AA(rankA+1:n,n+rankC+1:2*n) = MATMUL( Q3H, AA(rankA+1:n,n+rankC+1:2*n) )
   TEMP(1:n-rankA,1:n-rankC) = AA(rankA+1:n,n+rankC+1:two_n)
   CALL ZGEMM( 'NoTrans', 'NoTrans', n-rankA, n-rankC, n-rankA,                &
               CONE, Q3H, n-rankA, TEMP, two_n,                                &
               CZERO, AA(rankA+1,n+rankC+1), two_n )

   ! AA(1:rankA,1:n+rankC) = AA(1:rankA,1:n+rankC)*P3;
   ! AA(1:rankA,1:n+rankC) = MATMUL( AA(1:rankA,1:n+rankC), P3 )
   TEMP(1:rankA,1:n+rankC) = AA(1:rankA,1:n+rankC)
   CALL ZGEMM( 'NoTrans', 'NoTrans', rankA, n+rankC, n+rankC,                  &
               CONE, TEMP, two_n, P3, n+rankC, CZERO, AA(1,1), two_n )

   ! AA(n+1:n+rankC,1:n+rankC) = AA(n+1:n+rankC,1:n+rankC)*P3;
   ! AA(n+1:n+rankC,1:n+rankC) = MATMUL( AA(n+1:n+rankC,1:n+rankC), P3 )
   TEMP(1:rankC,1:n+rankC) = AA(n+1:n+rankC,1:n+rankC)
   CALL ZGEMM( 'NoTrans', 'NoTrans', rankC, n+rankC, n+rankC,                  &
               CONE, TEMP, two_n, P3, n+rankC, CZERO, AA(n+1,1), two_n )

   ! BB(:,1:n+rankC) = BB(:,1:n+rankC)*P3;
   ! BB(:,1:n+rankC) = MATMUL( BB(:,1:n+rankC), P3 )
   TEMP(:,1:n+rankC) = BB(:,1:n+rankC)
   CALL ZGEMM( 'NoTrans', 'NoTrans', two_n, n+rankC, n+rankC,                  &
               CONE, TEMP, two_n, P3, n+rankC, CZERO, BB(1,1), two_n )

   DEALLOCATE( TEMP, STAT = iwarn_STAT )

   ALLOCATE( ip(rankC+n-rankA), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! p = [(n+1:n+rankC) (rankA+1:n)];
   ip(1:rankC+n-rankA) = (/ (j, j = n+1,n+rankC), (j, j = rankA+1,n) /)
   AA(rankA+1:n+rankC,:) = AA(ip,:)
   BB(rankA+1:n+rankC,:) = BB(ip,:)

   DEALLOCATE( ip, STAT = iwarn_STAT )

   ALLOCATE( ip(n+rankC), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! p = [(n-rankA+1:n+rankC) (1:n-rankA)];
   ip(1:n+rankC) = (/ (j, j = n-rankA+1,n+rankC), (j, j = 1,n-rankA) /)
   AA(1:two_n,1:n+rankC) = AA(:,ip)
   BB(1:two_n,1:n+rankC) = BB(:,ip)

   DEALLOCATE( ip, STAT = iwarn_STAT )

   ALLOCATE( V(two_n,two_n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! V = [Z3 zeros(n+rankC,n-rankC); zeros(n-rankC,n+rankC) eye(n-rankC)];

   ! Form T = ( I  0  )*( P3  0 ).  P3 is  n+rankC by n+rankC  and
   !          ( 0  QC ) (  0  I )   QC is  n by n.
   V(1:two_n,1:two_n) = CZERO
   V(1:n+rankC,1:n+rankC) = P3
   DO j = n+rankC+1, two_n
      v(j,j) = CONE
   END DO

   DEALLOCATE( P3, STAT = iwarn_STAT )

   lwork = -1
   CALL ZUNMQR( 'Left', 'NoConjTrans', n, n, n, C, ldC, tauC, V(n+1,1), two_n, &
                dummy, lwork, INFO )
   lwork = dummy(1)
   ALLOCATE( work(lwork), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! V(n+1:two_n,:) = QC*V(n+1:two_n,:);
   CALL ZUNMQR( 'Left', 'NoConjTrans', n, n, n, C, ldC, tauC, V(n+1,1), two_n, &
                work, lwork, INFO )

   DEALLOCATE( work, STAT = iwarn_STAT )

   ! Now form V = T*( 0  I  0 )
   !                ( I  0  0 )
   !                ( 0  0  I )

   ! ind1 = (1:n+rankC);
   ! ind2 = [(n-rankA+1:n+rankC) (1:n-rankA)];
   ! V(:,ind1) = V(:,ind2);
   V(:,(/ (j, j = 1,n+rankC) /)) = V(:,(/ (j, j = n-rankA+1,n+rankC),          &
                                          (j, j = 1,n-rankA) /))
   IF( lvec )THEN
      ALLOCATE( Q(two_n,two_n), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! Form QC^H in Q(1:n:1:n)
      Q(1:n,1:n) = CZERO
      DO j = 1, n
         Q(j,j) = CONE
      END DO

      IF( nC /= ZERO )THEN
         lwork = -1
         CALL ZUNMQR( 'Left', 'ConjTrans', n, n, n, C, ldC, tauC, Q, two_n,    &
                      dummy, lwork, INFO )
         lwork = dummy(1)
         ALLOCATE( work(lwork), STAT = iwarn_STAT )
         IF( iwarn_STAT /= 0 )THEN
            INFO = -999
            iwarn = iwarn_STAT
            RETURN
         END IF

         CALL ZUNMQR( 'Left', 'ConjTrans', n, n, n, C, ldC, tauC, Q, two_n,    &
                      work, lwork, INFO )

         DEALLOCATE( work, STAT = iwarn_STAT )
      END IF

      ! Form the orthonormal matrix  Q = ( I   0    0  0 )*( QA^H   0   )
      !                                  ( 0   0    I  0 ) (  0    QC^H )
      !                                  ( 0  Q3^H  0  0 )
      !                                  ( 0   0    0  I )

      ! Q = [[QA(:,1:rankA) QA(:,rankA+1:n)*Q3]^H zeros(n);
      !       zeros(n)   QC'];
      ! ind1 = (rankA+1:n+rankC);
      ! ind2 = [(n+1:n+rankC) (rankA+1:n)];
      ! Q(ind1,:) = Q(ind2,:);
      ! Put QC^H in first to prevent overwriting that part of Q
      Q(rankA+1:rankA+rankC,n+1:two_n) = Q(1:rankC,1:n)
      Q(rankC+n+1:two_n,n+1:two_n) = Q(rankC+1:n,1:n)
      Q(1:rankA,1:n) = QAH(1:rankA,1:n)
      Q(1:rankA,n+1:two_n) = CZERO
      Q(rankA+1:rankA+rankC,1:n) = CZERO

      ! Q(rankA+rankC+1:rankC+n,1:n) = MATMUL( Q3H, QAH(rankA+1:n,1:n) )
      CALL ZGEMM( 'NoTrans', 'NoTrans', n-rankA, n, n-rankA,                   &
                  CONE, Q3H, n-rankA, QAH(rankA+1,1), n,                       &
                  CZERO, Q(rankA+rankC+1,1), two_n )

      Q(rankA+rankC+1:rankC+n,n+1:two_n) = CZERO
      Q(rankC+n+1:two_n,1:n) = CZERO
   END IF

   DEALLOCATE( Q3H, STAT = iwarn_STAT )
END IF

DEALLOCATE( jpvtA, STAT = iwarn_STAT )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Eigenvalue/eigenvector computation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! (1:iqz) is the active part of the linearization
iqz = rankC+rankA

! Assign the zero and infinite eigenvalues arising from rank deficient
! A and/or C.  The first iqz eigenvalues are computed later, by a call to ZGGEV,
! or ZGGES
IF( .NOT.rev )THEN
   alpha(iqz+1:n+rankC) = ONE
   beta(iqz+1:n+rankC) = ZERO
   alpha(n+rankC+1:two_n) = ZERO
   beta(n+rankC+1:two_n) = ONE
ELSE
   alpha(iqz+1:n+rankC) = ZERO
   beta(iqz+1:n+rankC) = ONE
   alpha(n+rankC+1:two_n) = ONE
   beta(n+rankC+1:two_n) = ZERO
END IF

IF( rvec )THEN
   ldVSR = two_n
ELSE
   ldVSR = 1
END IF
IF( lvec )THEN
   ldVSL = two_n
ELSE
   ldVSL = 1
END IF
ALLOCATE( VSL(ldVSL,ldVSL), VSR(ldVSR,ldVSR), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = -999
   iwarn = iwarn_STAT
   RETURN
END IF

! For the non-deflated eigenvalues use the QZ algorithm
! on active part of the big matrix pair (AA,BB).
ALLOCATE( rwork(MAX(1,8*iqz) ), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = -999
   iwarn = iwarn_STAT
   RETURN
END IF

IF( (iqz == two_n).OR.(.NOT.lvec) )THEN
   lwork = -1
   CALL ZGGEV( jobVL, jobVR, iqz, AA, two_n, BB, two_n, alpha, beta,           &
               VSL, ldVSL, VSR, ldVSR, dummy, lwork, rdummy, INFO )
   lwork = dummy(1)
   ALLOCATE( work(lwork), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   CALL ZGGEV( jobVL, jobVR, iqz, AA, two_n, BB, two_n, alpha, beta,           &
               VSL, ldVSL, VSR, ldVSR, work, lwork, rwork, INFO )
   IF( INFO /= 0 )THEN
      iwarn = INFO
      INFO = 3
      RETURN
   END IF

   DEALLOCATE( AA, BB, work, STAT = iwarn_STAT )
ELSE
   ! Left eigenvectors required and iqz < 2*n
   lwork = -1
   CALL ZGGES( 'Vec', jobVR, 'NoSort', ZSELCT, iqz, AA, two_n, BB, two_n,      &
               sdim, alpha, beta, VSL, ldVSL, VSR, ldVSR, dummy, lwork, rwork, &
               ldummy, INFO )
   lwork = dummy(1)
   lwork = MAX( lwork, 2*iqz )
   ALLOCATE( work(lwork), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! Compute the Schur form of the non-deflated part of the pencil, no
   ! eigenvectors at this stage.  Left and right Schur vectors are computed
   ! though
   CALL ZGGES( 'Vec', jobVR, 'NoSort', ZSELCT, iqz, AA, two_n, BB, two_n,      &
               sdim, alpha, beta, VSL, ldVSL, VSR, ldVSR, work, lwork, rwork,  &
               ldummy, INFO )
   IF( INFO /= 0 )THEN
      iwarn = INFO
      INFO = 2
      RETURN
   END IF

   IF( rvec )THEN
      ! Compute the right eigenvectors
      CALL ZTGEVC( 'Right', 'BackTrans', ldummy, iqz, AA, two_n, BB, two_n,    &
                   VSL, ldVSL, VSR, ldVSR, iqz, k, work, rwork, INFO )

   END IF
END IF
DEALLOCATE( rwork, work, STAT = iwarn_STAT )

IF( lvec.OR.rvec )THEN
   ALLOCATE( absa(two_n), alphan(two_n), rbeta(two_n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! Normalize the first iqz eigenvalues (alpha, beta) to lie on the
   ! unit circle
   CALL ZLANAB( iqz, alpha, beta, alphan, rbeta, absa, INFO )
   IF( INFO > 0 )THEN
      iwarn = INFO
      INFO = -999
   END IF

   IF( iqz < two_n )THEN
      rbeta(iqz+1:two_n) = beta(iqz+1:two_n)
   END IF
END IF

! Right eigenvector recovery
IF( rvec )THEN
   IF( ((scal == 1 .OR. scal == 2).AND.(tau < toltau)).OR.(rankC < n) )THEN
      ! The first vector always gives a good backward error

      ! VR(1:n,1:iqz) = MATMUL( V(1:n,1:iqz), VSR(1:iqz,1:iqz) )
      CALL ZGEMM( 'NoTrans', 'NoTrans', n, iqz, iqz, CONE, V, two_n,           &
                  VSR, ldVSR, CZERO, VR, ldVR )

      ! Normalize the vectors so that norm2(y) = 1 and the first
      ! element of largest absolute value is real and positive
      CALL ZLASGE( 2, n, iqz, VR, ldVR )

   ELSE
      ! ( scal == 0 .OR. scal == 3 .OR. scal == 4 .OR. tau >= toltau ).AND.
      ! (rankC == n)
      ! which also implies that  rankA = n

      ! Need to choose the vector with smallest backward error; compute
      ! first set of vectors
      ALLOCATE( TEMP(n,two_n), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! VSR(1:n,1:2*n) = MATMUL( V(1:n,1:2*n), VSR(1:2*n,1:2*n) )
      CALL ZGEMM( 'NoTrans', 'NoTrans', n, two_n, two_n, CONE, V, two_n,       &
                  VSR, ldVSR, CZERO, TEMP, n )

      VSR(1:n,1:two_n) = TEMP

      DEALLOCATE( TEMP, V, STAT = iwarn_STAT )

      lwork = -1
      CALL ZUNMQR( 'Left', 'ConjTrans', n, two_n, n, C, ldC, tauC,             &
                   VSR(n+1,1), ldVSR, dummy, lwork, INFO )
      lwork = dummy(1)
      ALLOCATE( work(lwork), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! Form QC^H*VSR(n+1:2*n,1:2*n) (with the result in VSR(n+1:2*n,1:2*n))
      CALL ZUNMQR( 'Left', 'ConjTrans', n, two_n, n, C, ldC, tauC,             &
                   VSR(n+1,1), ldVSR, work, lwork, INFO )

      ! Solve RC*X = VSR(n+1:2*n,1:2*n) for X
      ! (returning X in VSR(n+1:2*n,1:2*n))
      CALL ZTRSM( 'Left', 'Upper', 'NoTrans', 'NonUnit', n, two_n,             &
                  CONE, C, ldC, VSR(n+1,1), ldVSR )

      ! Now compute the second set of vectors
      ! VSR(n+1:2*n,1:2*n) = MATMUL(PC, VSR(n+1:2*n,1:2*n) )
      VSR(n+1:two_n,1:two_n) = VSR(n+jpvtC(:),1:two_n)

      DEALLOCATE( work, STAT = iwarn_STAT )

      ! Now choose the vectors with smallest backward error.  The vectors
      ! are returned in VR, are normalized and the element of largest
      ! absolute value is real and positive
      CALL ZLAG3V( n, two_n, alphan, rbeta, Acopy, n, B, ldB, Ccopy, n,        &
                   nA, nB, nC, VSR, ldVSR, VSR(n+1,1), ldVSR, VR, ldVR, INFO )
      IF( INFO > 0 )THEN
         iwarn = INFO
         INFO = -999
         RETURN
      END IF

      DEALLOCATE( VSR, STAT = iwarn_STAT )
   END IF

   ! Right eigenvectors corresponding to deflated eigenvalues.  Firstly, those
   ! coming from rank deficient C.
   IF( rankC < n )THEN
      IF( rankC == 0 )THEN
         VR(1:n,n+1:two_n) = CZERO
         DO j = 1, n
            VR(j,n+j) = ONE
         END DO
      ELSE
         lwork = -1
         CALL ZTZRZF( rankC, n, C, ldC, tauC, dummy, lwork, INFO )
         CALL ZUNMRZ( 'Right', 'ConjTrans', n, n, rankC, n-rankC, C, ldC,      &
                      tauC, PC, n, dummy(2), lwork, INFO )
         lwork = MAX( REAL( dummy(1) ), REAL( dummy(2) ) )
         ALLOCATE( work(lwork), STAT = iwarn_STAT )
         IF( iwarn_STAT /= 0 )THEN
            INFO = -999
            iwarn = iwarn_STAT
            RETURN
         END IF

         ! Reduce the upper trapezoidal matrix RC to upper triangular form
         CALL ZTZRZF( rankC, n, C, ldC, tauC, work, lwork, INFO )

         ! Transform PC, so that PC contains the (deflated) eigenvectors
         CALL ZUNMRZ( 'Right', 'ConjTrans', n, n, rankC, n-rankC, C, ldC,      &
                      tauC, PC, n, work, lwork, INFO )

         DEALLOCATE( work, STAT = iwarn_STAT )

         VR(1:n,n+rankC+1:two_n) = PC(:,rankC+1:n)
      END IF

      ! The vectors are already normalized, make the first element of largest
      ! absolute value positive
      CALL ZLASGE( 1, n, n-rankC, VR(1,n+rankC+1), ldVR)

   END IF
   ! Then those coming from rank deficient A.
   IF( rankA < n )THEN
      IF( rankA == 0 )THEN
         VR(1:n,rankC+1:n+rankC) = CZERO
         DO j = 1, n
            VR(j,rankC+j) = ONE
         END DO
      ELSE
         lwork = -1
         CALL ZTZRZF( rankA, n, A, ldA, tauA, dummy, lwork, INFO )
         CALL ZUNMRZ( 'Right', 'ConjTrans', n, n, rankA, n-rankA, A, ldA,      &
                      tauA, PA, n, dummy(2), lwork, INFO )
         lwork = MAX( REAL( dummy(1) ), REAL( dummy(2) ) )
         ALLOCATE( work(lwork), STAT = iwarn_STAT )
         IF( iwarn_STAT /= 0 )THEN
            INFO = -999
            iwarn = iwarn_STAT
            RETURN
         END IF

         ! Reduce the upper trapezoidal matrix RA to upper triangular form
         CALL ZTZRZF( rankA, n, A, ldA, tauA, work, lwork, INFO )

         ! Transform PA, so that PA contains the (deflated) eigenvectors
         CALL ZUNMRZ( 'Right', 'ConjTrans', n, n, rankA, n-rankA, A, ldA,      &
                      tauA, PA, n, work, lwork, INFO )

         DEALLOCATE( work, STAT = iwarn_STAT )

         VR(1:n,rankA+rankC+1:n+rankC) = PA(:,rankA+1:n)
      END IF

      ! The vectors are already normalized, make the first element of largest
      ! absolute value positive
      CALL ZLASGE( 1, n, n-rankA, VR(1,rankA+rankC+1), ldVR )

   END IF
END IF
DEALLOCATE( jpvtC, PA, PC, tauA, tauC, STAT = iwarn_STAT )

! Left eigenvectors recovery
IF( lvec )THEN
   ALLOCATE( TEMP(two_n,two_n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   ! If iqz = 2*n we have already computed the eigenvectors of the linearization
   IF( iqz < two_n )THEN
      ALLOCATE( rwork(4*n), select(two_n), work(4*n), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! VSL(1:iqz,1:iqz) = TRANSPOSE( CONJG( VSL(1:iqz,1:iqz) ) )
      ! AA(1:iqz,iqz+1:2*n) = MATMUL( VSL(1:iqz,1:iqz), AA(1:iqz,iqz+1:2*n) )
      ! BB(1:iqz,iqz+1:2*n) = MATMUL( VSL(1:iqz,1:iqz), BB(1:iqz,iqz+1:2*n) )
      ! Q(1:iqz,:) = MATMUL( VSL(1:iqz,1:iqz), Q(1:iqz,:) )
      TEMP(1:iqz,1:two_n-iqz) = AA(1:iqz,iqz+1:two_n)
      CALL ZGEMM( 'ConjTrans', 'NoTrans', iqz, two_n-iqz, iqz,                 &
                  CONE, VSL, ldVSL, TEMP, two_n, CZERO, AA(1,iqz+1), two_n )
      TEMP(1:iqz,1:two_n-iqz) = BB(1:iqz,iqz+1:two_n)
      CALL ZGEMM( 'ConjTrans', 'NoTrans', iqz, two_n-iqz, iqz,                 &
                  CONE, VSL, ldVSL, TEMP, two_n, CZERO, BB(1,iqz+1), two_n )
      TEMP(1:iqz,1:two_n) = Q(1:iqz,:)
      CALL ZGEMM( 'ConjTrans', 'NoTrans', iqz, two_n, iqz, CONE, VSL, ldVSL,   &
                  TEMP, two_n, CZERO, Q, two_n )

      select(1:iqz) = .TRUE.
      select(iqz+1:two_n) = .FALSE.
      VSL(1:two_n,1:iqz) = CZERO

      CALL ZTGEVC( 'Left', 'Selected', select, two_n, AA, two_n, BB, two_n,    &
                   VSL, ldVSL, dummy, 1, iqz, j, work, rwork, INFO )

      DEALLOCATE( AA, BB, rwork, select, work, STAT = iwarn_STAT )
   END IF

   ! VSL(1:2*n,1:iqz) = MATMUL( TRANSPOSE( CONJG( Q ) ), VSL(1:2*n,1:iqz) )
   TEMP(1:two_n,1:iqz) = VSL(:,1:iqz)
   CALL ZGEMM( 'ConjTrans', 'NoTrans', two_n, iqz, two_n, CONE, Q, two_n,      &
               TEMP, two_n, CZERO, VSL, two_n )

   DEALLOCATE( TEMP, STAT = iwarn_STAT )

   ! We now have two potential left eigenvectors in each column of VSL.
   ! Also, normalize the vectors so that norm2(y) = 1 and the first
   ! element of largest absolute value is real and positive
   IF( ( scal == 1 ).OR.( scal == 2 ).AND.( tau < toltau ) )THEN
      DO j = 1, iqz
         IF( ABS( alphan(j) ) < rbeta(j) )THEN
            VL(1:n,j) = VSL(n+1:two_n,j)
         ELSE
            VL(1:n,j) = VSL(1:n,j)
         END IF
      END DO

      CALL ZLASGE( 2, n, iqz, VL, ldVL )

   ELSE
      ! scal /= 1 or 2, or tau >= toltau

      ! Now choose the vectors with smallest backward error.  The vectors
      ! are returned in VL, are normalized and the element of largest
      ! absolute value is real and positive
      CALL ZLAG3V( n, iqz, alphan, rbeta, Acopy, n, B, ldB, Ccopy, n,        &
                   nA, nB, nC, VSL, ldVSL, VSL(n+1,1), ldVSL, VL, ldVL, INFO )
      IF( INFO > 0 )THEN
         iwarn = INFO
         INFO = -999
         RETURN
      END IF

   END IF

   ! Compute left eigenvectors of deflated zero and infinite eigenvalues.
   IF( rankA < n )THEN
      VL(1:n,iqz+1:n+rankC) = TRANSPOSE( CONJG( QAH(rankA+1:n,:) ) )

      ! The vectors are already normalized, make the element of largest
      ! absolute value real and positive
      CALL ZLASGE( 1, n, n+rankC-iqz, VL(1,iqz+1), ldVL )

   END IF
   IF( rankC < n )THEN
      VL(1:n,n+rankC+1:two_n) = QC

      ! The vectors are already normalized, make the element of largest
      ! absolute value positive
      CALL ZLASGE( 1, n, n-rankC, VL(1,n+rankC+1), ldVL )

      DEALLOCATE( QC, STAT = iwarn_STAT )
   END IF

END IF
DEALLOCATE( QAH, VSL, STAT = iwarn_STAT )

IF( lvec )THEN
   DEALLOCATE( Q, STAT = iwarn_STAT )
END IF

IF( rev )THEN
   ! Since we solved the reverse problem, need to reciprocate the eigenvalues
   ! (and keep beta real and non-negative),  The deflated eigenvalues were
   ! dealt with earlier
   DO j = 1, iqz
      IF( alpha(j) == CZERO )THEN
         alpha(j) = beta(j)
         beta(j) = ZERO
      ELSE
         tmp = ABS( alpha(j) )
         alpha(j) = ( CONJG( alpha(j) )/tmp )*beta(j)
         beta(j) = tmp
      END IF
   END DO
END IF

! If scaling took place, scale non-deflated eigenvalues back.
IF( scal > 0 )THEN
   alpha(1:iqz) = gg*alpha(1:iqz)
END IF

IF( lvec.OR.rvec )THEN
   IF( .not.rev )THEN
      A(1:n,1:n) = Acopy
      C(1:n,1:n) = Ccopy
   ELSE
      A(1:n,1:n) = Ccopy
      C(1:n,1:n) = Acopy
   END IF

   DEALLOCATE( Acopy, Ccopy, STAT = iwarn_STAT )

   IF( (gg /= ONE).OR.(delta /= ONE) )THEN
      ! Scale back A and C
      A(1:n,1:n) = (ONE/sa)*A(1:n,1:n)
      C(1:n,1:n) = (ONE/delta)*C(1:n,1:n)
   END IF

   ! Recover nA, nB and nC
   nA = g2
   nB = g1
   nC = g0
END IF

IF( (gg /= ONE).OR.(delta /= ONE) )THEN
   ! Scale back B
   B(1:n,1:n) = (ONE/sb)*B(1:n,1:n)
END IF

IF( sense > 0 )THEN
   ! Normalize the eigenvalues (alpha, beta) to lie on the unit circle.
   ! Note that the eigenvalues have been scaled, so this is a different
   ! normalization to earlier
   CALL ZLANAB( iqz, alpha, beta, alphan, rbeta, absa, INFO )
   IF( INFO > 0 )THEN
      iwarn = INFO
      INFO = -999
   END IF

   IF( .NOT.rev )THEN
      absa(iqz+1:n+rankC) = ONE
      absa(n+rankC+1:two_n) = ZERO
   ELSE
      absa(iqz+1:n+rankC) = ZERO
      absa(n+rankC+1:two_n) = ONE
   END IF
END IF

! Eigenvalue condition numbers
IF( ( sense == 1 ).OR.( sense > 4 ) )THEN
   ALLOCATE( abnu(3,two_n), DAX(n,iqz), TEMP(n,iqz), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   abnu(1,1:two_n) = ( absa**2 )*nA
   abnu(2,1:two_n) = absa*rbeta*nB
   abnu(3,1:two_n) = ( rbeta**2 )*nC

!  The first iqz (= rankA + rankC) condition numbers
   IF( nA /= ZERO )THEN
      DO j = 1, iqz
         TEMP(1:n,j) = ( TWO*alphan(j)*rbeta(j) )*VR(1:n,j)
      END DO

      CALL ZGEMM( 'NoTrans', 'NoTrans', n, iqz, n, CONE, A, ldA, TEMP, n,      &
                  CZERO, DAX, n )

   ELSE
      DAX(1:n,1:iqz) = CZERO
   END IF

   IF( nB /= ZERO )THEN
      DO j = 1, iqz
         TEMP(1:n,j) = ( (rbeta(j) + absa(j))*(rbeta(j) - absa(j)) )*VR(1:n,j)
      END DO

      CALL ZGEMM( 'NoTrans', 'NoTrans', n, iqz, n, CONE, B, ldB, TEMP, n,      &
                  CONE, DAX, n )

   END IF

   IF( nC /= ZERO )THEN
      DO j = 1, iqz
         TEMP(1:n,j) = ( -TWO*CONJG( alphan(j) )*rbeta(j) )*VR(1:n,j)
      END DO

      CALL ZGEMM( 'NoTrans', 'NoTrans', n, iqz, n, CONE, C, ldC, TEMP, n,      &
                  CONE, DAX, n )

   END IF

   DO j = 1, iqz
      asc = DOT_PRODUCT( VL(1:n,j), DAX(:,j) )
      IF( asc == CZERO )THEN
         s(j) = inf
      ELSE
         s(j) = DNRM2( 3, abnu(1,j), 1 )/ABS( asc )
      END IF
   END DO

   DEALLOCATE( DAX, TEMP, STAT = iwarn_STAT )

   IF( iqz < two_n )THEN
      ALLOCATE( DAX(n,two_n-iqz), STAT = iwarn_STAT )
      IF( iwarn_STAT /= 0 )THEN
         INFO = -999
         iwarn = iwarn_STAT
         RETURN
      END IF

      ! Condition numbers corresponding to the deflated eigenvalues,
      ! for which either alpha = 0, or beta = 0, so that
      ! DAX = B*VR(1:n,iqz+1:2*n)

      CALL ZGEMM( 'NoTrans', 'NoTrans', n, two_n-iqz, n, CONE, B, ldB,         &
                  VR(1,iqz+1), ldVR, CZERO, DAX, n )

      DO j = iqz+1, two_n
         asc = DOT_PRODUCT( VL(1:n,j), DAX(:,j-iqz) )
         IF( asc == CZERO )THEN
            s(j) = inf
         ELSE
            s(j) = DNRM2( 3, abnu(1,j), 1 )/ABS( asc )
         END IF
      END DO

      DEALLOCATE( DAX, STAT = iwarn_STAT )
   END IF

   DEALLOCATE( abnu, STAT = iwarn_STAT )
END IF

IF( sense > 1 )THEN
   ALLOCATE( p(two_n), STAT = iwarn_STAT )
   IF( iwarn_STAT /= 0 )THEN
      INFO = -999
      iwarn = iwarn_STAT
      RETURN
   END IF

   tmp = loctol*MAX( nA, nB, nC )
   p(1:two_n) = nA*absa**2 + nB*absa*rbeta + nC*rbeta**2
   IF( ( nA <= tmp ).OR.( nC <= tmp ) )THEN
      DO j = 1, two_n
         IF( p(j) <= tmp ) THEN
            p(j) = ONE
         END IF
      END DO
   END IF
END IF

IF( sense > 0 )THEN
   DEALLOCATE( absa, STAT = iwarn_STAT )
END IF

! Backward error of right eigenpairs
IF( (sense == 3).OR.(sense == 4).OR.(sense == 6).OR.(sense == 7) )THEN

   CALL ZLAG3B( 'Right', n, rev, rankA, rankC, alphan, rbeta, p,               &
                A, ldA, B, ldB, C, ldC, nA, nB, nC, VR, ldVR, beVR, INFO )

END IF

! Backward error of left eigenpairs
IF( (sense == 2).OR.(sense == 4).OR.(sense == 5).OR.(sense == 7) )THEN

   CALL ZLAG3B( 'Left', n, rev, rankA, rankC, alphan, rbeta, p,                &
                A, ldA, B, ldB, C, ldC, nA, nB, nC, VL, ldVL, beVL, INFO )

END IF

IF( scal > 1 )THEN
   DEALLOCATE( p, STAT = iwarn_STAT )
END IF

! Scale back to the original A, B, C
IF( lvec.OR.rvec )THEN
   DEALLOCATE( alphan, rbeta, STAT = iwarn_STAT )

   IF( d /= ONE )THEN
      A(1:n,1:n) = d*A(1:n,1:n)
      C(1:n,1:n) = d*C(1:n,1:n)
   END IF
END IF
IF( d /= ONE )THEN
   B(1:n,1:n) = d*B(1:n,1:n)
END IF

RETURN

END SUBROUTINE ZG3EVX
