!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Thomas R. Cameron* and Nikolas I. Steckley#
!                                                               
!   *Dept. Mathematics and Computer Science, Davidson College
!   #Dept. DiscoverOrg LLC.
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Last modified 21 December 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test Program for random initial approximations
!   
!   The convex hull is computed is Andrew's algorithm: monotone
!	chain algorithm.
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Inpute Parameter
!   
!   1) Start Degree, default 100
!   
!   2) End Degree, default 6400
!   
!   3) Iterations, default 10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM testStart
IMPLICIT NONE
! Input Variables
INTEGER :: startDegree, endDegree, iter, flag
CHARACTER(len=32) :: arg

! Compute Variables	
INTEGER :: clock, clock_rate, clock_start, clock_stop, deg, it
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: alpha, p, er, ei
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: timeMethods

! External Subroutines
EXTERNAL :: daruv, dsstart

! Set start and end degree, and iterations
CALL GET_COMMAND_ARGUMENT(1,arg,status=flag)
IF(flag==0) THEN
    READ(arg, '(I10)') startDegree
ELSE
	startDegree=100
ENDIF
CALL GET_COMMAND_ARGUMENT(2,arg,status=flag)
IF(flag==0) THEN
    READ(arg, '(I10)') endDegree
ELSE
    endDegree=6400
ENDIF
CALL GET_COMMAND_ARGUMENT(3,arg,status=flag)
IF(flag==0) THEN
    READ(arg, '(I10)') iter
ELSE
    iter=10
ENDIF

! Open Results File
OPEN(UNIT=1,FILE="results/testStart.csv")
WRITE(1,'(A)') 'Degree, Cameron Start Time, Bini Start Time'
! Allocate Elapsed Time storage arrays
ALLOCATE(timeMethods(iter,2))
! Main Loop
deg=startDegree
DO WHILE(deg<=endDegree)
    WRITE(1,'(I10)',advance='no') deg
    WRITE(1,'(A)',advance='no') ','
	
	DO it=1,iter
		! allocate polynomial and eigenvalue storage
		ALLOCATE(p(deg+1), alpha(deg+1), er(deg), ei(deg))
		! create random poly and its moduli
		CALL daruv(deg+1,p)
		alpha = dabs(p)
		! call Cameron start
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        CALL SYSTEM_CLOCK(count=clock_start)
		CALL start(alpha, deg, er, ei)
        CALL SYSTEM_CLOCK(count=clock_stop)
        timeMethods(it,1)=DBLE(clock_stop-clock_start)/DBLE(clock_rate)
		! call Bini start
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        CALL SYSTEM_CLOCK(count=clock_start)
		CALL dsstart(alpha, deg, er, ei)
        CALL SYSTEM_CLOCK(count=clock_stop)
        timeMethods(it,2)=DBLE(clock_stop-clock_start)/DBLE(clock_rate)
		! deallocate
		DEALLOCATE(p, alpha, er, ei)
	ENDDO
	
    WRITE(1,'(ES15.2)',advance='no') sum(timeMethods(:,1))/iter
    WRITE(1,'(A)',advance='no') ','
    WRITE(1,'(ES15.2)',advance='no') sum(timeMethods(:,2))/iter
    WRITE(1,'(A)')
	! update degree
	deg = 2*deg
ENDDO
! Deallocate Elapsed Time storage arrays
DEALLOCATE(timeMethods)
! Close Results File
CLOSE(1)

 CALL EXECUTE_COMMAND_LINE('python testStart.py')
CONTAINS

SUBROUTINE start(alpha, deg, er, ei)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: alpha(*)
DOUBLE PRECISION, INTENT(OUT)       :: er(*), ei(*)
!local scalars
INTEGER                             :: c, i, j, k, nzeros
DOUBLE PRECISION                    :: ang, r, th
!local arrays
INTEGER, DIMENSION(deg+1)           :: h
DOUBLE PRECISION, DIMENSION(deg+1)  :: a
!intrinsic procedures
INTRINSIC                           :: dcos, dexp, dlog, dsin
!parameters
DOUBLE PRECISION, PARAMETER         :: zero=0.0D0
DOUBLE PRECISION, PARAMETER         :: pi2 = 6.283185307179586D0, sigma = 0.7D0

! compute log(alpha)
DO i=1,deg+1
	IF(alpha(i)>zero) THEN
		a(i)=dlog(alpha(i))
	ELSE
		a(i)=-1.0D30
	ENDIF
ENDDO
! compute upper convex hull
CALL conv_hull(deg+1,a,h,c)
! compute initial estimates
k=0; th=pi2/deg
DO i=c-1,1,-1
	nzeros = h(i)-h(i+1)
	r = (alpha(h(i+1))/alpha(h(i)))**(1.0D0/nzeros)
	ang = pi2/nzeros
	DO j=1,nzeros
		er(k+j) = r*dcos(ang*j+th*h(i)+sigma)
		ei(k+j) = r*dsin(ang*j+th*h(i)+sigma)
	ENDDO
	k = k+nzeros
ENDDO
END SUBROUTINE start

SUBROUTINE conv_hull(n, a, h, c)
IMPLICIT NONE
! scalar argument
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(INOUT) :: c
! array argument
INTEGER, INTENT(INOUT)	:: h(*)
DOUBLE PRECISION, INTENT(IN) :: a(*)
! local scalars
INTEGER :: i
!parameters
DOUBLE PRECISION, PARAMETER :: zero=0.0D0

! build upper envelop of convex hull
 c=0
DO i=n,1,-1
	DO WHILE(c>=2 .and. cross(h,a,c,i)<=zero)
		c = c-1
	END DO
	c = c+1
	h(c)=i
END DO
END SUBROUTINE conv_hull

DOUBLE PRECISION FUNCTION cross(h, a, c, i)
IMPLICIT NONE
! scalar arguments
INTEGER, INTENT(IN) :: c, i
! array arguments
INTEGER, INTENT(IN) :: h(*)
DOUBLE PRECISION, INTENT(IN) :: a(*)

cross = (a(i)-a(h(c-1)))*(h(c)-h(c-1))-(a(h(c))-a(h(c-1)))*(i-h(c-1))
RETURN
END FUNCTION cross
	
END PROGRAM testStart
