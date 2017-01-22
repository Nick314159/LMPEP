PROGRAM tri_test_driver
USE environment
USE dgtlmpep_subroutines
USE qep3dlag
IMPLICIT NONE

!=======VARIABLES=======
!Common
INTEGER :: clock, clock_rate, clock_start, clock_stop, i, j, m, n, k, info
INTEGER :: degree, startDegree, maxDegree, jumpFactor
INTEGER :: mxit, iter, itermx, imax
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: timeStats, z
CHARACTER(LEN=100) :: arg
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, pdl, pd, pdu, xr, xi, yr, yi
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlangt
EXTERNAL :: dlangt
!------------------

!get information
CALL GETARG(1,arg)
READ(arg, *) startDegree
CALL GETARG(2,arg)
READ(arg, *) maxDegree
CALL GETARG(3, arg)
READ (arg,'(I10)') n
CALL GETARG(4, arg)
READ (arg,'(I10)') m
!CALL GETARG(3,BUFFER)
!READ(BUFFER, *) jumpFactor
jumpFactor=2

ALLOCATE(timeStats(m,2))
	
!Create iseed
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
iseed(4)=iseed(4)+1
ENDIF

OPEN(UNIT=1,FILE=resultsDir//"outputTri.csv")
WRITE(1, '(A)',  advance='no') 'DEGREE,  '
WRITE(1, '(A)',  advance='no') 'DGTLMPEP TIME,'
WRITE(1, '(A)',  advance='no') 'QEP3D TIME'
WRITE(1, *)
degree = startDegree
DO WHILE (degree < maxDegree)
  WRITE(1, '(i6)', advance='no') degree
  WRITE(1, '(A)', advance='no') ', '
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  DO j=1,m
    !create tridiagnol case
    ALLOCATE(work(n+n*(degree+1)), x(n), p(n,n*(degree+1)))
    DO i=1,degree+1
      CALL dlarnv(2, iseed, n, x)
      CALL dlagge(n, n, 1, 1, x, p(1,n*(i-1)+1), n, iseed, work, info)
    ENDDO
    DEALLOCATE(work, x)
    
    !store tridiagonal structure
    ALLOCATE(pdl(n-1,degree+1), pd(n,degree+1), pdu(n-1,degree+1))
    DO k=1,degree+1
      !pdl
      DO i=1,n-1
        pdl(i,k)=p(i+1,(k-1)*n+i)
      ENDDO
      !pd
      DO i=1,n
        pd(i,k)=p(i,(k-1)*n+i)
      ENDDO
      !pdu
      DO i=2,n
        pdu(i-1,k)=p(i-1,(k-1)*n+i)
      ENDDO
    ENDDO
    DEALLOCATE(p)
     
    
! Recursive subroutine for the evaluation of the eigenvalues of a 
! tridiagonal hyperbolic quadratic eigenvalue problem 
! Q(x) = x**2 tridiag(a2,a,a2) + x*tridiag(b2,b,b2) + tridiag(c2,c,c2)
!
! --------- output    ------
! z : eigenvalues
! iterall  : all iterations
! iterlast : iterations in the last conquer phase
! itermax  : maximum number of iterations in recursive calls
!
! --------- input     ------
! a,b,c      : diagonal of A, B, C, respectively
! a2, b2, c2 : codiagonal of A, B, C, respectively
! n          : matrix dimension
! maxit : maximum possible number of iterations
! mode  : method for the Newton correction: 1 (3-term), 2 (QR), 3 (LU)	

    ! maximal number of iteration
    mxit = 400	
    iter = 0
    itermx = 500
    !solve QEP3D
    ALLOCATE(z(2*n), zcx(2*n))
    CALL reigenl(a,au,b,bu,c,cu,n,z,mxit,iter,itermx,imax)
    DEALLOCATE(z,zcx)
    
    !solve dgtlmpep
    ALLOCATE(xr(n,n*degree), xi(n,n*degree), yr(n,n*degree), yi(n,n*degree))
    ALLOCATE(berr(n*degree), er(n*degree), ei(n*degree), ncoeff(degree+1), cond(n*degree), ferr(n*degree))
    DO i=1,degree+1
      ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
    ENDDO
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL dgtlm(pdl, pd, pdu, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, degree, n)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
  
    timeStats(j,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
 
    !deallocate
    DEALLOCATE(xr, xi, yr, yi)
    DEALLOCATE(berr, cond, ferr, er, ei, ncoeff, pdl, pd, pdu)
    
    !compute backward error and condition number
    CALL dposterrcond(pdl, pd, pdu, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, degree, n)
  
  END DO

  !=======SAVE RESULTS=======

  !dgtlmpep ----------
  WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,2))/m
  WRITE(1, '(A)', advance='no') ', '
  !------------------

  !qep3d -----------
  WRITE(1, '(20G15.4)', advance='no') SUM(timeStats(:,1))/m
  !------------------
  WRITE(1, *)
  degree = jumpFactor * degree
END DO
CLOSE(UNIT=1)
END PROGRAM tri_test_driver
