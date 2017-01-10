PROGRAM accuracy_driver
USE environment
USE dgelmpep_subroutines
IMPLICIT NONE

INTEGER ::  clock, clock_rate, clock_start, clock_stop, d, i, n, k, info, iwarn
REAL(dp) :: norm
CHARACTER (LEN=64), DIMENSION(27) :: tests
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, xr, xi, yr, yi
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, x
REAL(dp), DIMENSION(:), ALLOCATABLE :: alphar, alphai, beta, s, beVl, beVR
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: VL, VR
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlange, dznrm2
EXTERNAL :: dlange, dznrm2


tests(1) = 'bicycle.txt'
tests(2) = 'bilby.txt'
tests(3) = 'butterfly.txt'
tests(4) = 'cd_player.txt'
tests(5) = 'closed_loop.txt'
tests(6) = 'damped_beam.txt'
tests(7) = 'dirac.txt'
tests(8) = 'gen_hyper2.txt'
tests(9) = 'gen_tantipal2.txt'
tests(10) = 'gen_tpal2.txt'
tests(11) = 'hospital.txt'
tests(12) = 'intersection.txt'
tests(13) = 'metal_strip.txt'
tests(14) = 'mirror.txt'
tests(15) = 'mobile_manipulator.txt'
tests(16) = 'omnicam1.txt'
tests(17) = 'omnicam2.txt'
tests(18) = 'relative_pose_5pt.txt'
tests(19) = 'relative_pose_6pt.txt'
tests(20) = 'SKIP'!shaft.txt'
tests(21) = 'sleeper.txt'
tests(22) = 'speaker_box.txt'
tests(23) = 'spring_dashpot.txt'
tests(24) = 'spring.txt'
tests(25) = 'wing.txt'
tests(26) = 'wiresaw1.txt'
tests(27) = 'wiresaw2.txt'

CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF

OPEN(UNIT=1,FILE=resultsDir//"outputAccuracy.csv")
WRITE(1, '(A)',  advance='no') 'Problem, '
WRITE(1, '(A)',  advance='no') 'Max LM-BERR, '
WRITE(1, '(A)',  advance='no') 'Max LM-FERR, '
WRITE(1, '(A)',  advance='no') 'Max QM-BERR, '
WRITE(1, '(A)',  advance='no') 'Max QM-FERR, '
WRITE(1, *)

DO k=1,27
  !open problem
  IF(tests(k) .NE. 'SKIP') THEN
    OPEN(UNIT=2,FILE=problemsDir//tests(k))
    PRINT*, 'Testing '//tests(k)//'...'
    WRITE(1, '(A)', advance='no') tests(k)
    WRITE(1, '(A)', advance='no') ', '
    !read in size and degree from file
    READ(2,*) n
    READ(2,*) d
    !read in scalar polynomial
    ALLOCATE(p(n,n*(d+1)))
    READ(2,*) p
    CLOSE(UNIT=2)
  ELSE
    WRITE(*, '(A, I2)') 'Skipping test ',k
    GOTO 10
  ENDIF
  !solve problem using Laguerre's Method
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), er(n*d), ei(n*d), cond(n*d), ferr(n*d))
  ALLOCATE(ncoeff(d+1))
  DO i=1,d+1
    ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
  ENDDO
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NR')
  CALL SYSTEM_CLOCK(COUNT=clock_stop) 
  !bacward error, condition number for Laguerre's Method
  CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
  !Write results
  WRITE(1,'(ES15.2)', advance='no') MAXVAL(berr)
  WRITE(1, '(A)', advance='no') ', '
  WRITE(1,'(ES15.2)', advance='no') MAXVAL(ferr)
  WRITE(1, '(A)', advance='no') ', '
  !solve problem using QUADEIG
  IF(d==2) THEN
    ALLOCATE(alphar(2*n),alphai(2*n),beta(2*n))
    ALLOCATE(VL(n,2*n),VR(n,2*n),s(2*n),beVL(2*n),beVR(2*n))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL DG3EVX(1, 'V', 'V', 0, 0,                                    &
                   n, p(1,2*n+1), n, p(1,n+1), n, p(1,1), n,          &
                   alphar, alphai, beta, VL, n, VR, n, s, beVL, beVR, &
                   iwarn, info)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    !backward error, condition number for QUADEIG
    i=1
    DO WHILE(i<=2*n)
      IF(DABS(alphai(i))<eps) THEN
        !eigenvalues
        er(i)=alphar(i)/beta(i); ei(i)=zero
        !eigenvectors
        xr(:,i)=VR(:,i); xi(:,i)=zero
        yr(:,i)=VL(:,i); yi(:,i)=zero
        !normalize
        norm=dznrm2(n,DCMPLX(xr(:,i)),1)
        xr(:,i)=xr(:,i)/norm
        norm=dznrm2(n,DCMPLX(yr(:,i)),1)
        yr(:,i)=yr(:,i)/norm
        i=i+1
      ELSE
        !eigenvalues
        er(i)=alphar(i)/beta(i); ei(i)=alphai(i)/beta(i)
        er(i+1)=alphar(i)/beta(i); ei(i+1)=-alphai(i)/beta(i)
        !eigenvectors
        xr(:,i)=VR(:,i); xi(:,i)=VR(:,i+1)
        yr(:,i)=VL(:,i); yi(:,i)=VL(:,i+1)
        xr(:,i+1)=VR(:,i); xi(:,i+1)=-VR(:,i+1)
        yr(:,i+1)=VL(:,i); yi(:,i+1)=-VL(:,i+1)
        !normalize
        norm=dznrm2(n,DCMPLX(xr(:,i),xi(:,i)),1)
        xr(:,i)=xr(:,i)/norm; xi(:,i)=xi(:,i)/norm
        xr(:,i+1)=xr(:,i+1)/norm; xi(:,i+1)=xi(:,i+1)/norm
        norm=dznrm2(n,DCMPLX(yr(:,i),yi(:,i)),1)
        yr(:,i)=yr(:,i)/norm; yi(:,i)=yi(:,i)/norm
        yr(:,i+1)=yr(:,i+1)/norm; yi(:,i+1)=yi(:,i+1)/norm
        i=i+2
      ENDIF
    ENDDO
    CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr,d,n)
    !Write results
    WRITE(1,'(ES15.2)', advance='no') MAXVAL(berr)
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1,'(ES15.2)', advance='no') MAXVAL(ferr)
    WRITE(1, '(A)', advance='no') ', '
    !Deallocate
    DEALLOCATE(alphar,alphai,beta,VL,VR,s,beVL,beVR)
  ELSE
    WRITE(1, '(A)', ADVANCE='no') 'NA, NA, '
  ENDIF
  !Deallocate
  DEALLOCATE(p, xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)
  WRITE(1,*)
  10 CONTINUE
ENDDO

END PROGRAM accuracy_driver
