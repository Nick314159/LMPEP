! Demo program for the computation of the eigenvalues of the 
! tridiagonal quadratic eigenvalue problems using the Ehrlich-Aberth
! and the Laguerre method
! 
! Reference: B. Plestenjak, Numerical methods for the tridiagonal
! hyperbolic quadratic eigenvalue problem, Preprint 967, IMFM, 
! Ljubljana, 2005
!
! Bor Plestenjak
! bor.plestenjak@fmf.uni-lj.si
! Department of Mathematics
! University of Ljubljana
! November 2005

PROGRAM main

	USE qep3dea
	USE qep3deacx
	USE qep3dlag

	IMPLICIT NONE

	INTEGER								:: mode, n, i, neg, detsgn
	INTEGER, PARAMETER					:: ATTEMPTS = 200
	REAL(dp), ALLOCATABLE, DIMENSION(:)	:: a, b, c, au, bu, cu, al, bl, cl, d
	REAL(dp), ALLOCATABLE, DIMENSION(:)	:: z
	COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: zcx
	REAL(4)								:: T1
	INTEGER								:: mxit, iter, itermx, imax
	
!	Data is read from standard input
!	First parameter is mode
!	- 1 : symmetric hyperbolic tridiagonal E-A + 3-term
!	- 2 : symmetric hyperbolic tridiagonal E-A + QR
!	- 3 : symmetric hyperbolic tridiagonal E-A + LUMV
!	- 4 : symmetric hyperbolic tridiagonal Laguerre + 3-term
!	- 5 : tridiagonal E-A + 3-term
!	- 6 : tridiagonal E-A + QR
!	- 7 : tridiagonal E-A + LUMV
!   Second parameter is the size n
!   Then we have vectors a, a2, b, b2, c, c2 (mode 1-4)
!   or vectors a, au, al, b, bu, bl, c, cu, cl 
!	in column forms
!	
!	Check the prepared data examples data_n_j
!	
!	Output is on standard output in the form
!	time
!	iter, itemmx, imax
!	real(z(1)) imag(z(1))
!	...
!	real(z(2n)) imag(z(2n))

	READ(*,*) mode
	READ(*,*) n
    ALLOCATE(a(n), b(n), c(n), au(n), al(n), bu(n), bl(n), cu(n), cl(n))
    ALLOCATE(z(2*n), zcx(2*n))
    IF (mode<5) THEN
		READ(*,*) ( a(i), i=1,n)
		READ(*,*) ( au(i), i=1,n-1)
		READ(*,*) ( b(i), i=1,n)
		READ(*,*) ( bu(i), i=1,n-1)
		READ(*,*) ( c(i), i=1,n)
		READ(*,*) ( cu(i), i=1,n-1)
	ELSE
		READ(*,*) ( a(i), i=1,n)
		READ(*,*) ( au(i), i=1,n-1)
		READ(*,*) ( al(i), i=1,n-1)
		READ(*,*) ( b(i), i=1,n)
		READ(*,*) ( bu(i), i=1,n-1)
		READ(*,*) ( bl(i), i=1,n-1)
		READ(*,*) ( c(i), i=1,n)
		READ(*,*) ( cu(i), i=1,n-1)
		READ(*,*) ( cl(i), i=1,n-1)
	END IF
	
	z = zero
	zcx = zero

	
	mxit = 400	! maximal number of iteration
	iter = 0
	itermx = 500
	T1 = SECNDS(0.0)	
	! Real Ehrlich-Aberth
	IF (mode>=1 .AND. mode<=3) THEN
		CALL reigen(a,au,b,bu,c,cu,n,z,mxit,iter,itermx,imax,mode)
		WRITE(*,*) SECNDS(T1)
		WRITE(*,*) iter, itermx, imax
		WRITE(*,'(E25.15E3," ",E25.15E3)') (z(i), 0.d0, i=1,2*n)
	END IF

	! Laguerre
	IF (mode==4) THEN
		CALL reigenl(a,au,b,bu,c,cu,n,z,mxit,iter,itermx,imax)
		WRITE(*,*) SECNDS(T1)
		WRITE(*,*) iter, itermx, imax
		WRITE(*,'(E25.15E3," ",E25.15E3)') (z(i), 0.d0, i=1,2*n)
	END IF

	! Complex Ehrlich Aberth
	IF (mode>=5 .AND. mode<=7) THEN
		CALL reigencx(a,au,al,b,bu,bl,c,cu,cl,n,zcx,mxit,iter,itermx,imax,mode-4)
		WRITE(*,*) SECNDS(T1)
		WRITE(*,*) iter, itermx, imax
		WRITE(*,'(EN25.15E3," ",EN25.15E3)') (REAL(zcx(i)), AIMAG(zcx(i)), i=1,2*n)
	END IF

END PROGRAM main




