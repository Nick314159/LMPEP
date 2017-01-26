PROGRAM tri_test_driver2
USE environment
USE dgtlmpep_subroutines
USE eigensolve
IMPLCIIT NONE

!===Variables===
!LMPEP

!EIGENSOLVE
INTEGER                  :: n,i,j,m
REAL(dp),DIMENSION(:)    :: a,s,rad ,ss,tt,cond
COMPLEX(dp),DIMENSION(:) :: z 
ALLOCATABLE              :: a, s, z, rad,ss,tt,cond
REAL(dp)                 :: iter ,aaa,bbb,theta,h
CHARACTER(len=20)        :: filename
REAL                     :: ru

!size
n=100
!allocate
ALLOCATE(a(n), s(n), z(n), cond(n), rad(n))
DO j=1,8
  IF(j==1)THEN
!!! TEST 1
    DO i=1,n
      a(i)=i*(-1.d0)**(i/8)
      s(i)=((-1.d0)**i)/i
    END DO
  END IF
  IF(j==2)THEN
!!! TEST 2
    DO i=1,n
      a(i)= 10*(-1.d0)**(i/8)
      s(i)=i*(-1)**(i/9)
    END DO
  END IF
  IF(j==3)THEN
!!! TEST 3
    DO i=1,n
      a(i)=i
      s(i)=(n-i+1)
    END DO
  END IF
  IF(j==4)THEN
!!! TEST 4
    DO i=1,n
      a(i)=1.d0*(-1)**i
      s(i)=20.d0*(-1)**(i/5)
    END DO
  END IF
  IF(j==5)THEN
!!! TEST 5. 
    DO i=1,n
      a(i)=(-1)**(i/4)*10.d0**(5*(-1)**i)
      s(i)=(-1)**(i/3)
    END DO
  END IF
  IF(j==6)THEN
!!! TEST 6
    DO i=1,n
      a(i)=2
      s(i)=1
    END DO
  END IF
  IF(j==7)THEN
!!! TEST 7
    DO i=1,n
      a(i)=1.d0/i+1.d0/(n-i+1)
      s(i)=(1.d0/i)*(-1)**(i/9)
    END DO
  END IF
  IF(j==8)THEN
!!! TEST 8
    DO i=1,n
      a(i)=i*(-1)**(i/13)*(-1)**(i/5)
      s(i)=(-1)**(i/11)*(n-i+1)**2.d0
    END DO
  END IF
!!! Compute eigenvalues using eigen
  CALL eigen(n,a,s,z,con)
ENDDO


END PROGRAM tri_test_driver2
