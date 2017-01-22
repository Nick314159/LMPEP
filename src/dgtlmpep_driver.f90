PROGRAM dgtlmpep_driver
USE dgtlmpep_subroutines
IMPLICIT NONE

!local scalars
INTEGER :: c, clock, clock_rate, clock_start, clock_stop, d, i, info, j, n
CHARACTER (LEN=32) f
!local arrays
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x, y, tau, co, si
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, pdl, pd, pdu, xr, xi, yr, yi
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlangt
EXTERNAL :: dlangt
!file location (where problem files are stored)
CHARACTER(*), PARAMETER :: fileplace1="/home/thomas/Documents/FORTRAN/LMPEP2/PROBLEMS/REAL/"
!file location (where results are stored)
CHARACTER(*), PARAMETER :: fileplace2="/home/thomas/Documents/FORTRAN/LMPEP2/RESULTS/"

!create iseed, used in dlarnv and dlagge
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

!create or open problem
IF(c==1) THEN
  !read in degree
  PRINT*, 'input size and degree'
  READ*, n
  READ*, d
  ALLOCATE(work(n+n*(d+1)), x(n), p(n,n*(d+1)))
  DO i=1,d+1
    CALL dlarnv(2, iseed, n, x)
    CALL dlagge(n, n, 1, 1, x, p(1,n*(i-1)+1), n, iseed, work, info)
  ENDDO
  DEALLOCATE(work, x)
ELSEIF(c==2) THEN
  !read in file name
  PRINT*, 'input file name'
  READ*, f
  OPEN(UNIT=1,FILE=fileplace1//f)
  !read in size and degree from file
  READ(1,*) n
  READ(1,*) d
  !read in scalar polynomial
  ALLOCATE(p(n,n*(d+1)))
  READ(1,*) p
  CLOSE(UNIT=1)
ENDIF
!store tridiagonal structure
ALLOCATE(pdl(n-1,d+1), pd(n,d+1), pdu(n-1,d+1))
DO j=1,d+1
  !pdl
  DO i=1,n-1
    pdl(i,j)=p(i+1,(j-1)*n+i)
  ENDDO
  !pd
  DO i=1,n
    pd(i,j)=p(i,(j-1)*n+i)
  ENDDO
  !pdu
  DO i=2,n
    pdu(i-1,j)=p(i-1,(j-1)*n+i)
  ENDDO
ENDDO
DEALLOCATE(p)

!solve problem
ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
ALLOCATE(berr(n*d), er(n*d), ei(n*d), ncoeff(d+1), cond(n*d), ferr(n*d))
DO i=1,d+1
  ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
ENDDO
CALL SYSTEM_CLOCK(count_rate=clock_rate)
CALL SYSTEM_CLOCK(COUNT=clock_start)
CALL dgtlm(pdl, pd, pdu, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n)
CALL SYSTEM_CLOCK(COUNT=clock_stop)

!compute backward error and condition number
CALL dposterrcond(pdl, pd, pdu, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)

!save results
WRITE(*,'(20G15.4)') 'DGTLMPEP TIME =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)
WRITE(*,'(20G15.4)') 'MAX BERR      =', MAXVAL(berr)
WRITE(*,'(20G15.4)') 'MAX FERR      =', MAXVAL(ferr)

OPEN(UNIT=1,FILE=fileplace2//"EIGVALUES")
DO i=1,n*d
  WRITE(1,*) er(i), ei(i)
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=fileplace2//"BERR-COND-FERR")
DO i=1,n*d
  WRITE(1,*) berr(i), cond(i), ferr(i)
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=fileplace2//"R-EIGVECTORS")
DO j=1,n*d
  DO i=1,n
    WRITE(1,*) xr(i,j), xi(i,j)
  ENDDO
  WRITE(1,*) NEW_LINE('c')
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=fileplace2//"L-EIGVECTORS")
DO j=1,n*d
  DO i=1,n
    WRITE(1,*) yr(i,j), yi(i,j)
  ENDDO
  WRITE(1,*) NEW_LINE('c')
ENDDO
 CLOSE(UNIT=1)
!deallocate
DEALLOCATE(xr, xi, yr, yi)
DEALLOCATE(berr, cond, ferr, er, ei, ncoeff, pdl, pd, pdu)

END PROGRAM dgtlmpep_driver
