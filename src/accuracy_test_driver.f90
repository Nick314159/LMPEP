PROGRAM accuracy_driver
USE environment
USE dgeeam_subroutines
IMPLICIT NONE

!local scalars
INTEGER :: c, clock, clock_rate, clock_start, clock_stop, d, i, info, j, n, k
!local arrays
CHARACTER (LEN=64), DIMENSION(27) :: tests
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, xr, xi, yr, yi
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x
INTEGER, DIMENSION(4) :: iseed
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK

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
WRITE(1, '(A)',  advance='no') 'Problem,     '
WRITE(1, '(A)',  advance='no') 'Max BERR,     '
WRITE(1, '(A)',  advance='no') 'Max FERR,    '
WRITE(1, '(A)', advance='no') 'Time        '
WRITE(1, *)

DO k=1,27
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
  
    !solve problem using Laguerre's Method
    ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
    ALLOCATE(berr(n*d), er(n*d), ei(n*d), cond(n*d), ferr(n*d))
    ALLOCATE(ncoeff(d+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NR')
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
  
    !bacward error, condition number for Laguerre's Method
    CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
  
    !Write results
    WRITE(1,'(20G15.4)', advance='no') MAXVAL(berr)
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1,'(20G15.4)', advance='no') MAXVAL(ferr)
    WRITE(1,'(20G15.4)', advance='no') DBLE(clock_stop-clock_start)/DBLE(clock_rate)

    WRITE(1,*)

    DEALLOCATE(p, xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)
  ELSE
    WRITE(*, '(A, I2)') 'Skipping test ',k
  END IF
END DO
 CLOSE(UNIT=1)

END PROGRAM accuracy_driver
