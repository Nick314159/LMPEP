! Driver program for the normalized problem
! with refinement of the approximations in quadruple precision
PROGRAM driver
  USE eigensolve
  IMPLICIT NONE
  INTEGER                  :: n, i, j
  REAL(dp),DIMENSION(:)    :: a, s, cond 
  COMPLEX(dp),DIMENSION(:) :: z 
  real(qp),dimension(:)    :: qa, qs
  COMPLEX(qp),DIMENSION(:) :: qz 
  ALLOCATABLE              :: a, s, z, cond, qa, qs, qz
  REAL(dp)                 :: iter 
  CHARACTER(len=20)        :: filename
  REAL                     :: ru
  
  WRITE(*,*)"n="
  READ(*,*)n
  WRITE(*,*)"computing the eigenvalues of a normalized random problem..."
  
  ALLOCATE(a(n),s(n),z(n),cond(n), qa(n),qs(n),qz(n))
  WRITE(*,*)"computing the eigenvalues of a random problem..."
  ! Random matrix
  DO i=1,n
     CALL RANDOM_NUMBER(ru)
     a(i)=0.5-ru
     CALL RANDOM_NUMBER(ru)
     s(i)=0.5-ru
  END DO

  write(*,*)"computing eigenvalues..."
  CALL eigen(n,a,s,z,cond)
  write(*,*)"writing eigenvalues in the file 'eigenvalues'"
  OPEN(unit=2,file="eigenvalues")
  DO i=1,n
     WRITE(2,*)REAL(z(i)),AIMAG(z(i))
  END DO
  WRITE(*,*)"refining eigenvalues"
  qa=a
  qs=s
  qz=z
  call qaberth(n,qa,qs,qz)
   write(*,*)"writing refined eigenvalues in the file 'qeigenvalues'"
  OPEN(unit=12,file="qeigenvalues")
  DO i=1,n
     WRITE(12,*)REAL(qz(i)),AIMAG(qz(i))
  END DO
END PROGRAM driver
