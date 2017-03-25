PROGRAM dslmpep_driver
USE dslmpep_subroutines
IMPLICIT NONE

!local scalars
INTEGER :: c, clock, clock_rate, clock_start, clock_stop, d, i
CHARACTER (LEN=16) f
!local arrays
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, er, ei, p
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAXVAL, MOD, SYSTEM_CLOCK
!file location (where problem files are stored)
CHARACTER(*), PARAMETER :: fileplace1="/home/thomas/Documents/FORTRAN/LMPEP2/PROBLEMS/REAL/"
!file location (where results are stored)
CHARACTER(*), PARAMETER :: fileplace2="/home/thomas/Documents/FORTRAN/LMPEP2/RESULTS/"

!create iseed, used in dlarnv
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF

!random problem or stored problem?
PRINT*, 'input 1 for random problem, 2 for stored problem'
READ*, c

!create random or open problem
IF(c==1) THEN
  !read in degree
  PRINT*, 'input degree'
  READ*, d
  ALLOCATE(p(d+1))
  CALL dlarnv(2,iseed,d+1,p)
ELSEIF(c==2) THEN
  !read in file name
  PRINT*, 'input file name'
  READ*, f
  OPEN(UNIT=1,FILE=fileplace1//f)
  !read in degree from file
  READ(1,*) d
  !read in scalar polynomial
  ALLOCATE(p(d+1))
  READ(1,*) p
  CLOSE(UNIT=1)
ENDIF

!solve problem
ALLOCATE(berr(d),er(d),ei(d))
berr=zero
CALL SYSTEM_CLOCK(count_rate=clock_rate)
CALL SYSTEM_CLOCK(count=clock_start)
CALL dslm(p, er, ei, berr, d)
CALL SYSTEM_CLOCK(count=clock_stop)

!save results
WRITE(*,'(20G15.4)') 'DSLMPEP TIME =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)
WRITE(*,'(20G15.4)') 'AVG. BERR    =', SUM(berr)/d
WRITE(*,'(20G15.4)') 'MAX BERR     =', MAXVAL(berr)
OPEN(UNIT=1,FILE=fileplace2//"ROOTS-BERR")
DO i=1,d
  WRITE(1,*) er(i), ei(i), berr(i)
ENDDO
 CLOSE(UNIT=1)

!deallocate
DEALLOCATE(berr,er,ei,p)

END PROGRAM dslmpep_driver
