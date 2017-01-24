!! Driver file with the tests

PROGRAM driver
USE eigensolve
  IMPLICIT NONE
  INTEGER                  :: n,i,j,m
  REAL(dp),DIMENSION(:)    :: a,s,rad ,ss,tt,cond
  COMPLEX(dp),DIMENSION(:) :: z 
  ALLOCATABLE              :: a, s, z, rad,ss,tt,cond
  REAL(dp)                 :: iter ,aaa,bbb,theta,h
  CHARACTER(len=20)        :: filename
  REAL                     :: ru


  WRITE(*,*)"data from file: type 0, else type the test number (1-16)"
  write(*,*)"0- file"
  write(*,*)"1- eigenvalues along three lines"
  write(*,*)"2- real eigs plus eigenvalues close to the imag axis"
  write(*,*)"3- real eigs between 10^-4 and n"
  write(*,*)"4- eigs along 5 curves"
  write(*,*)"5- eigs in 5 clusters, of large and small moduli"
  write(*,*)"6- real eigs: tridiag [1,2,1]"
  write(*,*)"7- eigs along 16 lines from the origin"
  write(*,*)"8- eigs mostly real"
  write(*,*)"9- eigs with a fractal-like shape"
  write(*,*)"10- eigs in 5 tight clusters of large and small moduli"
  write(*,*)"11- Pasquini's: very ill condit eigs. n<20, a= ,b= "
  write(*,*)"12- eigs as in test 10 with exact input"
  write(*,*)"13- Liu matrix n=28: eig=0 multiplicity=n"
  write(*,*)"14- Liu matrix n=14: eig=0 multiplicity=n"
  write(*,*)"15- Clement matrix: integer eigenvalues "
  write(*,*)"16- Dorr matrix "
  write(*,*)">=17- Random "

  READ(*,*) j
  IF(j==0)THEN
     WRITE(*,*)"file name ="
     READ(*,*)filename
     OPEN(unit=2,file=filename)
     READ(2,*)n
     ALLOCATE(a(n), s(n))
     ALLOCATE(z(n),cond(n))
     DO i=1,n
        READ(2,*)a(i)
     END DO
     DO i=1,n
        READ(2,*)s(i)
     END DO
  ELSE
	if(j==13.or.j==14.or.j<0) goto 100
     WRITE(*,*)"type the dimension n="
     READ(*,*)n
     ALLOCATE(a(n), s(n))
     allocate(z(n),cond(n))
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
!!! * TEST 8
        DO i=1,n
           a(i)=i*(-1)**(i/13)*(-1)**(i/5)
           s(i)=(-1)**(i/11)*(n-i+1)**2.d0
        END DO
     END IF
     IF(j==9)THEN
!!! * TEST 9 The matrix {{-x,1,0},{1,x,1},{0,1,x}} is such that det =-x^3
!!!  the matrix obtained by concatenating two blocks equal to the previous
!!!  matrix has nonzero eigenvalues 
        DO i=1,n
           a(i)=1
           s(i)=1.d0
           IF(i>=n/2)s(i)=-1.d0
        END DO
     END IF
     IF(j==10)THEN
        DO i=1,n
           a(i)=-(10.d0)**(5*(-1)**i)
           s(i)=1
        END DO
        s(1)=-1
        a(n)=10.d0**5
     END IF
     IF(j==11)THEN
!!! Test 11:  Pasquini
        allocate(ss(n),tt(n))
        write(*,*)"a=";read(*,*)aaa
        write(*,*)"b=";read(*,*)bbb
        a(1)=-bbb/aaa
        ss(1)=-a(1)
        tt(1)=a(1)/(aaa+1)
        do i=2,n-1
           a(i)=-bbb*(aaa-2.d0)/((2*i+aaa-2.d0)*(2*i+aaa-4))
           ss(i)=bbb*(i+aaa-2)/((2*i+aaa-2)*(2*i+aaa-3))
           tt(i)=-bbb*i/((2*i+aaa-1.d0)*(2*i+aaa-2.d0))
        end do
        a(n)=-bbb*(aaa-2.d0)/((2*n+aaa-2.d0)*(2*n+aaa-4.d0))
        s=1.d0
        call normalize(n,a,ss,tt,s)
     END IF
!!! TEST 12
IF(j==12)THEN
        DO i=1,n
           a(i)= 2.d0**(20 *(-1)**(i+1))
           if(i<=n/2)then
              s(i)=1.d0
           else
              s(i)=-1.d0
           end if
        END DO
     END IF
!!! Test 13 Liu n=28
100  continue
     If(j==13)Then
        n=28
        allocate(a(n),s(n),z(n),cond(n))
        allocate(ss(n),tt(n))
        a=0
        a(7)=-1;a(8)=1;a(14)=-1;a(15)=1;a(21)=1;a(22)=-1
        s=1;ss=1;tt=1
        ss(1)=-1;ss(4)=-1;ss(6:8)=-1
        ss(10)=-1;ss(13:15)=-1;ss(18)=-1
        ss(20:22)=-1;ss(24)=-1;ss(27)=-1
        call normalize(n,a,ss,tt,s)
     end if
!!! Test 14 Liu n=14
     If(j==14)Then
        n=14
        allocate(a(n),s(n),z(n),cond(n))
        allocate(ss(n),tt(n))
        a=0
        a(7)=-1;a(8)=1;a(14)=0
        s=1;ss=1;tt=1
        ss(1)=-1;ss(4)=-1;ss(6:8)=-1
        ss(10)=-1;ss(13)=-1
        call normalize(n,a,ss,tt,s)
     end if
          If(j==-1)Then
        n=28
        allocate(a(n),s(n),z(n),cond(n))
        allocate(ss(n),tt(n))
        a=0
        a(7)=-1;a(8)=1;a(14)=-1;a(15)=1;a(21)=1;a(22)=-1
        s=1;ss=1;tt=1
        ss(1)=-1;ss(4)=-1;ss(6:8)=-1
        ss(10)=-1;ss(13:15)=-1;ss(18)=-1
        ss(20:22)=-1;ss(24)=-1;ss(27)=-1
        a(16)=64.d0
        call normalize(n,a,ss,tt,s)
     end if

!!! Test 15 Clement
     if(j==15)then
        a=0
        s=1
        allocate(ss(n),tt(n))
        tt=0
        do i=1,n-1
           tt(i)=i*(n-i)
           tt(i)=sqrt(tt(i))
        end do
        ss=tt
        call normalize(n,a,ss,tt,s)
     end if
     
!!! Test 16 Dorr
     If(j==16)then
        allocate(ss(n),tt(n))
        write(*,*)'theta='
        read(*,*)theta
        h=1.d0/(n+1)
        m=(n+1)/2
        theta=theta/h**2
        do i=1,m
           ss(i)=-theta
           tt(i)=ss(i)-(0.5d0-i*h)/h
           a(i)=-(ss(i)+tt(i))
        end do
        do i=m+1,n
           tt(i)=-theta
           ss(i)=tt(i)+(0.5d0-i*h)/h
           a(i)=-(ss(i)+tt(i))
        end do
        ss(1:n-1)=ss(2:n)
        s=1
        call normalize(n,a,ss,tt,s)
     end if
!!! TEST 17 Random
     IF(j>16)THEN
        DO i=1,n
           CALL RANDOM_NUMBER(ru)
           a(i)=0.5-ru
           CALL RANDOM_NUMBER(ru)
           s(i)=0.5-ru
        END DO
     END IF
 end if 

!!! Compute eigenvalues
  CALL eigen(n,a,s,z,cond)
  OPEN(unit=2,file="ploteigen")
  OPEN(unit=12,file="eig")
  OPEN(unit=3,file="cond")
  DO i=1,n
     WRITE(2,*)REAL(z(i)),AIMAG(z(i))
      WRITE(12,*)z(i)
    write(3,*)cond(i)
  END DO

!!! Remove the STOP if validation is required
  STOP
  ALLOCATE(rad(n))
  CALL validate(n,a,s,z,rad)
  OPEN(unit=3,file="radii")
  DO i=1,n
     WRITE(3,*)"z=",z(i),"rad=",10.d0**rad(i)
  END DO
END PROGRAM driver

