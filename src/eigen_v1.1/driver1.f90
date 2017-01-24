!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                       MODULE EIGENSOLVE                             !!!
!!!!!              by D.A. Bini, L. Gemignani, F. Tisseur                 !!!
!!!!!                              v. 1.1                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                Latest revision: March 2004                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample of driver program for the normalized problem
PROGRAM driver
  USE eigensolve
  IMPLICIT NONE
  INTEGER                  :: n, i, j
  REAL(dp),DIMENSION(:)    :: a, s, cond 
  COMPLEX(dp),DIMENSION(:) :: z 
  ALLOCATABLE              :: a, s, z, cond
  REAL(dp)                 :: iter 
  CHARACTER(len=20)        :: filename
  REAL                     :: ru
  
  WRITE(*,*)"n="
  READ(*,*)n
  WRITE(*,*)"computing the eigenvalues of a normalized random problem..."
  
  ALLOCATE(a(n),s(n),z(n),cond(n))
  WRITE(*,*)"computing the eigenvalues of a random problem..."
  ! Random matrix
  DO i=1,n
     CALL RANDOM_NUMBER(ru)
     a(i)=0.5-ru
     CALL RANDOM_NUMBER(ru)
     s(i)=0.5-ru
  END DO

  CALL eigen(n,a,s,z,cond)
  OPEN(unit=2,file="eigenvalues")
  DO i=1,n
     WRITE(2,*)REAL(z(i)),AIMAG(z(i))
  END DO
  WRITE(*,*)"...done. Output in the file 'eigenvalues'"
  OPEN(unit=12,file="cond")
  DO i=1,n
     WRITE(12,*)cond(i)
  END DO
  WRITE(*,*)"Condition number estimates written in the file 'cond'"
  
END PROGRAM driver
