! Driver program for the  unnormalized problem (T-zS)v=0
! a: diagonal entries of T
! b: superdiagonal entries of T
! c: subdiagonal entries of T
! s: diagonal entries of S
PROGRAM driver
USE eigensolve
  IMPLICIT NONE
  INTEGER                  :: n,i,j
  REAL(dp),DIMENSION(:)    :: a, s, rad, b, c, cond
  COMPLEX(dp),DIMENSION(:) :: z 
  ALLOCATABLE              :: a, s, z, rad, b, c, cond
  REAL(dp)                 :: iter 
  CHARACTER(len=20)        :: filename
  REAL                     :: ru

  WRITE(*,*)"n="
  READ(*,*)n
  ALLOCATE(a(n),b(n),c(n),s(n),z(n),cond(n))
  WRITE(*,*)"computing the eigenvalues of a random problem..."
  ! Random matrix
  DO i=1,n
     CALL RANDOM_NUMBER(ru)
     a(i)=0.5-ru
     CALL RANDOM_NUMBER(ru)
     s(i)=0.5-ru
     CALL RANDOM_NUMBER(ru)
     b(i)=0.5-ru
     CALL RANDOM_NUMBER(ru)
     c(i)=0.5-ru
  END DO
  WRITE(*,*)"Normalizing..."
  call normalize(n,a,b,c,s)
  WRITE(*,*)"...done. Now computing eigenvalues..."
  CALL eigen(n,a,s,z,cond)
  OPEN(unit=2,file="eigenvalues")
  DO i=1,n
     WRITE(2,*)z(i)
  END DO
  WRITE(*,*)"...done. Output in the file 'eigenvalues'"
  WRITE(*,*)"Validating the result..." 
  ALLOCATE(rad(n))
  CALL validate(n,a,s,z,rad)
  WRITE(*,*)"...done. Output in the file 'validated'"
  OPEN(unit=3,file="validated")
  DO i=1,n
     j=rad(i)+0.5
     WRITE(3,*)"z=",z(i),"rad=",10.d0**rad(i)
  END DO
END PROGRAM driver


