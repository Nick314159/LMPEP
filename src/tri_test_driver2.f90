PROGRAM tri_test_driver2
USE environment
USE dgtlmpep_subroutines2
USE eigensolve
IMPLICIT NONE

!===Variables===
!LMPEP
INTEGER :: clock, clock_rate, clock_start, clock_stop, d
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: er, ei, ncoeff, berr
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pd, pdl, pdu
!EIGENSOLVE
INTEGER                              :: n,i,j
REAL(dp),DIMENSION(:),ALLOCATABLE    :: a,s,cond
COMPLEX(dp),DIMENSION(:),ALLOCATABLE :: z
REAL                                 :: ru	
!external procedures
REAL(dp) :: dlangt
EXTERNAL :: dlangt

!create iseed, used in dgtlmpep
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF

!===Cost Complexity Tests===
PRINT*, 'Cost Complexity Tests'
n=10
d=1
DO WHILE(n<400)
PRINT*, 'size ', n
!!! Allocate
  ALLOCATE(a(n), s(n), z(n), cond(n*d))
  ALLOCATE(pdl(n-1,d+1), pd(n,d+1), pdu(n-1,d+1))
  ALLOCATE(ncoeff(d+1))
  ALLOCATE(berr(n*d), er(n*d), ei(n*d))
!!! Create random problem
  DO i=1,n
    CALL RANDOM_NUMBER(ru)
    a(i)=0.5-ru
    CALL RANDOM_NUMBER(ru)
    s(i)=0.5-ru
  ENDDO
!!! Compute eigenvalues using DGTLMPEP
  pdl(:,1)=one; pd(:,1)=a; pdu(:,1)=one
  pdl(:,2)=zero; pd(:,2)=-s; pdu(:,2)=zero
  !coefficient norm
  DO i=1,d+1
    ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
  ENDDO
  !solve problem
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL dgtlm(pdl,pd,pdu,er,ei,berr,ncoeff,iseed,d,n)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)
  PRINT*, 'DGTLMPEP time =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)

!!! Compute eigenvalues using EIGEN
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL eigen(n,a,s,z,cond)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)
  PRINT*, 'EIGEN time =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)

!!! Deallocate
  DEALLOCATE(a, s, z, cond)
  DEALLOCATE(pdl, pd, pdu)
  DEALLOCATE(ncoeff)
  DEALLOCATE(berr, er, ei)
!!! Update n
n=2*n
ENDDO

END PROGRAM tri_test_driver2
