! Module for the computation of the eigenvalues of the hyperbolic 
! tridiagonal real symmetric quadratic eigenvalue problem using the 
! Laguerre method
! 
! References: 
!
! [1] B. Plestenjak, Numerical methods for the tridiagonal
!     hyperbolic quadratic eigenvalue problem, Preprint 967, IMFM, 
!     Ljubljana, 2005
!
! [2] module EIGENSOLVE, by D.A. Bini, L. Gemignani, F. Tisseur for 
!     the generalized eigenvalue problem, V. 1.1, March 2004 
!
! Bor Plestenjak
! bor.plestenjak@fmf.uni-lj.si
! Department of Mathematics
! University of Ljubljana
! November 2005


MODULE qep3dlag

USE qepinit

IMPLICIT NONE

CONTAINS

!======================================================================
!	ddet3dsr3t
!======================================================================

SUBROUTINE simplesort(z,n,nmid)

! Reorders vector z of length n, where parts [1..nmid] and [nmid+1..n]
! are ordered

! --- arguments ---	
INTEGER							:: n
REAL(dp), DIMENSION(n)			:: z
! --- local ---
REAL(dp), DIMENSION(n)			:: t
INTEGER							:: inda, indb, i, nmid
REAL(dp)						:: tmp

inda = 1
indb = nmid+1

DO i=1,n
	
	IF ((inda<=nmid) .AND. (indb<=n)) THEN
		IF (z(inda)<z(indb)) THEN
			t(i) = z(inda)
			inda = inda+1
		ELSE
			t(i) = z(indb)
			indb = indb+1
		END IF
	ELSE 
		IF (inda>nmid) THEN
			t(i) = z(indb)
			indb = indb + 1
		ELSE 
			t(i) = z(inda)
			inda = inda + 1
		END IF
	END IF 

END DO

z = t

END subroutine simplesort

!======================================================================
!	ddet3dsr3t
!======================================================================

SUBROUTINE ddet3dsr3t(a,a2,b,b2,c,c2,n,x,d,df,ddf,neg,detsgn,mode)

! Computes the first and the second derivative of the determinant 
! of a real symmetric tridiagonal quadratic polynomial using 3-term 
! reccurences
!
! mode:
!		1 : neg + dsign (for bisection) and d
!		2 : df, ddf (for Laguerre's method) where d is known (input)
!		3 : df, ddf (for Laguerre's method) + neg + dsign 
!	 
! --- arguments ---	
INTEGER							:: n, neg, mode, detsgn
REAL(dp), DIMENSION(n)			:: a, b, c
REAL(dp), DIMENSION(n-1)		:: a2, b2, c2
REAL(dp)						:: x, df, ddf
! --- local ---
REAL(dp), DIMENSION(n)			:: t, dt, ddt, d, g, h
REAL(dp), DIMENSION(n-1)		:: t2, dt2, ddt2
INTEGER							:: i
    
IF (.NOT. mode==2) THEN

	df = zero
	ddf = zero
	detsgn = 0

	! The determinant of Q(x)

	t	= x**2*a  + x*b + c			! diagonal of Q(x)
	t2	= x**2*a2 + x*b2 + c2		! codiagonal of Q(x)

	d(1) = t(1)
	IF (d(1)==0) THEN
		d(1) = eps*(eps + ABS(x**2*a(1)) + ABS(x*b(1)) + ABS(c(1)))
	END IF
	d(2) = t(2)-t2(1)**2/d(1);
	IF (d(2)==0) THEN
		d(2) = eps*(eps + ABS(x**2*a(2)) + ABS(x*b(2)) + ABS(c(2)))
	END IF 
	DO i=2,n-1
		d(i+1) = t(i+1)-t2(i)**2/d(i)
		IF (d(i+1)==0) THEN 
			d(i+1) = eps*(eps + ABS(x**2*a(i+1)) + ABS(x*b(i+1)) + ABS(c(i+1))) 
		END IF
	END DO

	! Negative count and sign of the determinant
!	IF (mode==1) THEN
		neg = COUNT(d<0)
		IF (MOD(neg,2)==1) THEN 
			detsgn = -1
		ELSE
			detsgn = 1
		END IF
!	END IF
END IF

! The first and the second derivative of det(Q(x))
	
IF (mode>1) THEN

	! The first derivative
	
	dt	= 2*x*a + b					! diagonal of Q'(x)	
	dt2 = 2*x*a2 + b2				! codiagonal of Q'(x)
	g(1) = dt(1)/d(1) 
	g(2) = ( dt(2) + t(2)*g(1) - (2*t2(1)*dt2(1))/d(1) )/d(2)
	DO i=2,n-1
		g(i+1) = ( dt(i+1) + t(i+1)*g(i) &
			& - (2*t2(i)*dt2(i) + t2(i)**2*g(i-1))/d(i) )/d(i+1)
	END DO   
	df = g(n)

	! The second derivative

	ddt	= 2*a 				 		! diagonal of Q''(x)	
	ddt2 = 2*a2						! codiagonal of Q''(x)
	h(1) = ddt(1)/d(1)
	h(2) = ( ddt(2) + 2*dt(2)*g(1) + t(2)*h(1) &
			& -(2*dt2(1)**2 + 2*t2(1)*ddt2(1))/d(1) )/d(2)
	DO i=2,n-1
		h(i+1) = ( ddt(i+1) + 2*dt(i+1)*g(i) + t(i+1)*h(i) & 
			& -(2*dt2(i)**2 + 2*t2(i)*ddt2(i) &
			& + 4*t2(i)*dt2(i)*g(i-1) + t2(i)**2*h(i-1))/d(i) )/d(i+1)
	END DO
	ddf = h(n)

END IF

END SUBROUTINE ddet3dsr3t


!======================================================================
!	reigenl
!======================================================================

RECURSIVE SUBROUTINE reigenl(a, a2, b, b2, c, c2, n, z, maxit, & 
						&   iterall, iterlast, itermax)

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
!
! --------- arguments ------
INTEGER                    :: n, maxit, iterall, iterlast, itermax
REAL(dp), DIMENSION(n)     :: a, b, c
REAL(dp), DIMENSION(n-1)   :: a2, b2, c2
REAL(dp), DIMENSION(2*n)   :: z
!
! --------- local     ------
INTEGER                    :: i, nm, itermax1, itermax2, itermax3
INTEGER                    :: iter1, iter2, tmpiter
REAL(dp)                   :: discr, sqrtdiscr, pertb
REAL                       :: rnd

! for n=1 the eigenvalues are solution of a quadratic equation
    
IF (n==1) THEN
    discr = b(1)**2-4*a(1)*c(1)
    IF (discr>=0) THEN
		sqrtdiscr = SIGN(SQRT(discr),-b(1))
		z(1) = (- b(1) + sqrtdiscr )/(2*a(1))				 
		z(2) = c(1)/a(1)/z(1)
	ELSE
		! something is wrong, the eigenvalues should be real
		z(1) = -b(1)/(2*a(1));
		z(2) = -b(1)/(2*a(1));
	END IF
	iterall = 0
	iterlast = 0
	itermax = 0
    RETURN
 END IF
       
! problem of size n>1 is split into two subproblems

nm=(n+1)/2 
	   	
CALL reigenl(a(1:nm),a2(1:nm-1),b(1:nm),b2(1:nm-1),c(1:nm),c2(1:nm-1),& 
					& nm,z(1:2*nm),maxit,iter1,tmpiter,itermax)

CALL reigenl(a(nm+1:n),a2(nm+1:n-1),b(nm+1:n),b2(nm+1:n-1),c(nm+1:n),c2(nm+1:n-1),&
					& n-nm,z(2*nm+1:2*n),maxit,iter2,tmpiter,itermax)

! it is essential to sort the eigenvalues
CALL simplesort(z,2*n,2*nm)

! refine the approximations by the Laguerre iteration
CALL laguerre(a, a2, b, b2, c, c2, n, z, maxit, iterlast, itermax)
iterall = iter1 + iter2 + iterlast   

END SUBROUTINE reigenl

!======================================================================
!	laguerre
!======================================================================

SUBROUTINE laguerre(a,a2,b,b2,c,c2,n,z,maxit,iterall,iterlevel)

! This subroutine performs the Laguerre iterations for refining the 
! approximation z to the eigenvalue of
! Q(x) = x**2*tridiag(a2,a,a2) + x*tridiag(b2,b,b2) + tridiag(c2,c,c2)
!
! We use the following stopping criteria:
!	1) the correction is small
!	2) the inertia has changed
!
! --------- output    ------
! z : eigenvalues
! iterall : number of iterations required 
! iterlevel : maximum level of iteration
!
! --------- input     ------
! a,b,c      : diagonal of A, B, C, respectively
! a2, b2, c2 : codiagonal of A, B, C, respectively
! n          : matrix dimension
! z			 : initial approximations
! maxit : maximum possible number of iterations
!
! --------- arguments ------
INTEGER                     :: n, maxit, iterall, iterlevel, mode
REAL(dp), DIMENSION(n)      :: a, b, c
REAL(dp), DIMENSION(n-1)    :: a2, b2, c2
REAL(dp), DIMENSION(2*n)    :: z
!
! --------- local     ------
INTEGER                     :: i, j, m, nfound, razlog, lastchange
REAL(dp)                    :: corr, abcorr, newt, den, oldcorr 
LOGICAL,DIMENSION(2*n)      :: cont
REAL(dp), DIMENSION(2*n)    :: minres
INTEGER, DIMENSION(2*n)		:: kand
REAL(dp)					:: relchange, rnd, imen1, imen2
INTEGER						:: mrem
REAL(dp), DIMENSION(2*n)    :: newz
INTEGER						:: detsgn, neg, goal
REAL(dp)					:: x, df, ddf, xa, xb, xc, xd, razlika
REAL(dp), DIMENSION(n)		:: d
INTEGER						:: pm, iter, konec, oldpm, bisiter, xchg
LOGICAL						:: pogoj, pog1, pog2, pog3

m = 2*n
iterall = 0

DO i=1,m
	
	x = z(i)
	
	! the inertia is checked for the direction and corrected for the
	! secondary eigenvalues
	CALL ddet3dsr3t(a,a2,b,b2,c,c2,n,x,d,df,ddf,neg,detsgn,3)
	IF (i>n) THEN
		neg = m - neg
	END IF		
	
	! based on the inertia we have to go left (-1) or right (1) 
	IF (neg<i) THEN
		pm = 1
	ELSE
		pm = -1
	END IF

	! conditions for a good approximation
	! a) neg=i or neg=i-1
	! b) the sign of f'/f agrees with -pm
	! c) the sign of f' is (-1)^i
	pog1 = (ABS(neg-i)<=1)
	pog2 = (.NOT. df*pm>0)
	pog3 = ((-1)**(i)*df*detsgn >= 0)
	pogoj = (pog1 .AND. pog2 .AND. pog3)

	konec = 0

	! If the above conditions are not true, we do bisection
	bisiter = 0
	IF ((i>1) .AND. (i<m) .AND. (.NOT. i==n) .AND. (.NOT. i==n+1) .AND. ((.NOT. pogoj) .OR. (ABS(neg-i)>1))) THEN
		xa = MIN(x,z(i+pm))
		xb = MAX(x,z(i+pm))
		
		DO WHILE ( (.NOT. pogoj) .AND. (bisiter<5) .AND. (ABS(xa-xb)>eps*ABS(x)))
			bisiter =  bisiter + 1
			! first two steps are special, as we do not cut in half
			! the reason is that in many cases one of the ends of the
			! interval is very close to the eigenvalue we are looking for
			IF (bisiter==1) THEN
				x = (one-1.d-06)*xa + 1.d-06*xb
			ELSE IF (bisiter==2) THEN
				x = (one-1.d-06)*xb + 1.d-06*xa
			ELSE
				x = (xa + xb) / 2
			END IF

			CALL ddet3dsr3t(a,a2,b,b2,c,c2,n,x,d,df,ddf,neg,detsgn,3)
			IF (i>n) THEN
				neg = m - neg
			END IF		

	        IF (neg<i) THEN
				xa = x
			ELSE
				xb = x
			END IF
			IF (neg<i) THEN
				pm = 1
			ELSE
				pm = -1
			END IF

			pog1 = (ABS(neg-i)<=1)
			pog2 = (.NOT. df*pm>0)
			pog3 = ((-1)**(i)*df*detsgn >= 0)
			pogoj = (pog1 .AND. pog2 .AND. pog3)
			IF (ABS(xa-xb)<eps*ABS(x)) THEN
				konec = -1
			END IF

		END DO
	END IF

	! now we do Laguerre iteration
	iter = 0
	oldpm = pm

	DO WHILE ((konec==0) .AND. (iter<maxit)) 
		iter = iter + 1
		CALL ddet3dsr3t(a,a2,b,b2,c,c2,n,x,d,df,ddf,neg,detsgn,3)
		
		IF (i>n) THEN
			neg = m - neg
		END IF		
		
		IF (neg<i) THEN
			pm = 1
		ELSE
			pm = -1
		END IF
		pogoj = (df*pm<0)
		den = (m-1)*((m-1)*df**2-m*ddf)
		IF (den<0) THEN
			write(11,*) "WARNING, complex square root n=",n,"i=",i,"x=",x,"den=",den
			den = ABS(den)
		END IF 			
        imen1 = -df+SQRT(den);     
        imen2 = -df-SQRT(den);
        IF (((pm>0) .AND. (imen1>imen2)) .OR. ((pm<0) .AND. (imen1<imen2))) THEN
           corr = m/imen1
        ELSE
           corr = m/imen2
        END IF

		x = x + corr
		razlika = ABS(corr)/ABS(x)
		
		IF (razlika<eps) THEN
			konec = 1
		ELSE IF (.NOT. pm==oldpm) THEN
			konec = 2
		END IF
		
	END DO

	z(i) = x
	iterall = iterall + iter + bisiter
	IF (iter + bisiter >iterlevel) THEN
		iterlevel = iter
	END IF
END DO

END SUBROUTINE laguerre

END MODULE qep3dlag












