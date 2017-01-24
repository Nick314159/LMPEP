! Module for the computation of the eigenvalues of the 
! tridiagonal real quadratic eigenvalue problem using the 
! Ehrlich-Aberth method. Eigenvalue can be complex as well.
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

MODULE qep3deacx

USE util

IMPLICIT NONE

CONTAINS

!======================================================================
!	newtcorr3Tcx
!======================================================================

SUBROUTINE newtcorr3Tcx(a, au, al, b, bu, bl, c, cu, cl, n, x, df)

! Computes the Newton correction for the determinant of a real 
! tridiagonal quadratic polynomial matrix using
! 3-term reccurences
!
! f(x) = det(x^2*tridiag(al,a,au)+x*tridiag(bl,b,bu)+tridiag(cl,c,cu))
!
! --------- output    ------
! df : f'(x)/f(x)
!
! --------- input     ------
! a,b,c : diagonal of A, B, C, respectively
! al, bl, cl : lower codiagonal of A, B, C, respectively
! au, bu, cu : upper codiagonal of A, B, C, respectively
! n : matrix dimension
! x : complex value
!
! --------- arguments ------	
INTEGER						:: n
REAL(dp), DIMENSION(n)		:: a, b, c
REAL(dp), DIMENSION(n-1)	:: al, bl, cl, au, bu, cu
COMPLEX(dp)					:: x, df
!
! --------- local     ------
COMPLEX(dp), DIMENSION(n)	:: t, dt, d, g
COMPLEX(dp), DIMENSION(n-1)	:: tu, tl, dtu, dtl
INTEGER						:: i
    
! The determinant of Q(x)
	
t	= x**2*a  + x*b + c			! diagonal of Q(x)
tu	= x**2*au + x*bu + cu		! upper codiagonal of Q(x)
tl	= x**2*al + x*bl + cl		! lower codiagonal of Q(x)

d(1) = t(1)
IF (d(1)==0) THEN
	d(1) = eps*(eps + ABS(x**2*a(1)) + ABS(x*b(1)) + ABS(c(1)))
END IF
d(2) = t(2)-tu(1)*tl(1)/d(1);
IF (d(2)==0) THEN
	d(2) = eps*(eps + ABS(x**2*a(2)) + ABS(x*b(2)) + ABS(c(2)))
END IF 
DO i=2,n-1
	d(i+1) = t(i+1)-tu(i)*tl(i)/d(i)
	IF (d(i+1)==0) THEN 
		d(i+1) = eps*(eps + ABS(x**2*a(i+1)) + ABS(x*b(i+1)) + ABS(c(i+1))) 
	END IF
END DO

dt	= 2*x*a + b					! diagonal of Q'(x)	
dtu = 2*x*au + bu				! upper codiagonal of Q'(x)
dtl = 2*x*al + bl				! lower codiagonal of Q'(x)

g(1) = dt(1)/d(1) 
g(2) = (dt(2)+t(2)*g(1)-1/d(1)*(tu(1)*dtl(1)+tl(1)*dtu(1)))/d(2)
DO i=2,n-1
	g(i+1) = (dt(i+1)+t(i+1)*g(i)-(tu(i)*dtl(i) + tl(i)*dtu(i) + tu(i)*tl(i)*g(i-1))/d(i))/d(i+1)
END DO   
	
df = g(n)

END SUBROUTINE newtcorr3Tcx

!======================================================================
!	newtcorr3Tcx
!======================================================================

SUBROUTINE newtcorr3Tcx2(a, au, al, b, bu, bl, c, cu, cl, n, x, df)

! Computes the Newton correction for the determinant of a real 
! tridiagonal quadratic polynomial matrix using
! 3-term reccurences
!
! f(x) = det(x^2*tridiag(al,a,au)+x*tridiag(bl,b,bu)+tridiag(cl,c,cu))
!
! --------- output    ------
! df : f'(x)/f(x)
!
! --------- input     ------
! a,b,c : diagonal of A, B, C, respectively
! al, bl, cl : lower codiagonal of A, B, C, respectively
! au, bu, cu : upper codiagonal of A, B, C, respectively
! n : matrix dimension
! x : complex value
!
! --------- arguments ------	
INTEGER						:: n
REAL(dp), DIMENSION(n)		:: a, b, c
REAL(dp), DIMENSION(n-1)	:: al, bl, cl, au, bu, cu
COMPLEX(dp)					:: x, df
!
! --------- local     ------
COMPLEX(dp), DIMENSION(n)	:: t, dt, d, g
COMPLEX(dp), DIMENSION(n-1)	:: tu, tl, dtu, dtl
INTEGER						:: i
    
! The determinant of Q(x)
	
t	= x**2*a  + x*b + c			! diagonal of Q(x)
tu	= x**2*au + x*bu + cu		! upper codiagonal of Q(x)
tl	= x**2*al + x*bl + cl		! lower codiagonal of Q(x)

d(1) = t(1)
IF (d(1)==0) THEN
	d(1) = eps*(eps + ABS(x**2*a(1)) + ABS(x*b(1)) + ABS(c(1)))
END IF
d(2) = t(2)-tu(1)*tl(1)/d(1);
IF (d(2)==0) THEN
	d(2) = eps*(eps + ABS(x**2*a(2)) + ABS(x*b(2)) + ABS(c(2)))
END IF 
DO i=2,n-1
	d(i+1) = t(i+1)-tu(i)*tl(i)/d(i)
	IF (d(i+1)==0) THEN 
		d(i+1) = eps*(eps + ABS(x**2*a(i+1)) + ABS(x*b(i+1)) + ABS(c(i+1))) 
	END IF
END DO

dt	= 2*x*a + b					! diagonal of Q'(x)	
dtu = 2*x*au + bu				! upper codiagonal of Q'(x)
dtl = 2*x*al + bl				! lower codiagonal of Q'(x)

g(1) = dt(1)/d(1) 
g(2) = (dt(2)+t(2)*g(1)-1/d(1)*(tu(1)*dtl(1)+tl(1)*dtu(1)))/d(2)
DO i=2,n-1
	g(i+1) = (dt(i+1)+t(i+1)*g(i)-(tu(i)*dtl(i) + tl(i)*dtu(i) + tu(i)*tl(i)*g(i-1))/d(i))/d(i+1)
END DO   
	
df = g(n)

END SUBROUTINE newtcorr3Tcx2

!======================================================================
!	givenscx
!======================================================================

SUBROUTINE givenscx(x, y, c, s, r)

! Computes Givens rotation c,s,r such that
!
! computes c,s,r such that
!
! [c               s] [x] = [r]
! [-conj(s)  conj(c)] [y]   [0]
!
! with conj(c)*c + conj(s)*s = 1
!
! --------- arguments ------	
COMPLEX(dp)					:: x,y,c,s,r
!
! --------- local     ------
REAL(dp)					:: scale, absx, absy

absx = ABS(x)
absy = ABS(Y)
IF (absx .EQ. zero) THEN
    c = zero
	r = absy
    s = CONJG(y)/r
ELSE
    scale = absx + absy
    r = scale * SQRT((absx/scale)**2 + (absy/scale)**2)
	c = CONJG(x)/r;
	s = CONJG(y)/r
END IF

END SUBROUTINE givenscx
	
!======================================================================
!	newtcorrQRcx
!======================================================================

SUBROUTINE newtcorrQRcx(a, au, al, b, bu, bl, c, cu, cl, n, x, tr)

! computes derivative of the determinant of a symmetric tridiagonal 
! quadratic matrix polynomial using QR decomposition
!	 
! f(x) = det(x^2*tridiag(al,a,au)+x*tridiag(bl,b,bu)+tridiag(cl,c,cu))
!
! --------- output    ------
! tr : f'(x)/f(x)
!
! --------- input     ------
! a,b,c : diagonal of A, B, C, respectively
! al, bl, cl : lower codiagonal of A, B, C, respectively
! au, bu, cu : upper codiagonal of A, B, C, respectively
! n : matrix dimension
! x : complex value
!
! --------- arguments ------	
INTEGER						:: n
REAL(dp), DIMENSION(n)		:: a, b, c
REAL(dp), DIMENSION(n-1)	:: al, bl, cl, au, bu, cu
COMPLEX(dp)					:: x, tr
!
! --------- local     ------
COMPLEX(dp), DIMENSION(n)	:: alpha, bf, r, u, v, w
COMPLEX(dp), DIMENSION(n-1)	:: beta, gama, bg, bh, ps, s
COMPLEX(dp), DIMENSION(n-2)	:: t
COMPLEX(dp)					:: aa, g, phi
INTEGER						:: i
    
alpha = x**2*a  + x*b  + c		! diagonal of Q(x)
beta = x**2*al + x*bl + cl		! lower codiagonal of Q(x)
gama = x**2*au + x*bu + cu		! upper codiagonal of Q(x)
bf = 2*x*a   + b				! diagonal of Q'(x)	
bg = 2*x*al  + bl				! lower codiagonal of Q'(x)
bh = 2*x*au  + bu				! upper codiagonal of Q'(x)

aa = alpha(1)
g = gama(1)
u(1) = one
ps = zero

! Computes vectors r, s, t defining R and u and v defining Q.
DO i=1,n-1
	CALL givenscx(aa, beta(i), phi, ps(i), r(i))
	s(i) = -CONJG(ps(i))*(phi*g+ps(i)*alpha(i+1))
	aa = -CONJG(ps(i))*g + CONJG(phi)*alpha(i+1)
	v(i) = phi
	u(i+1) = CONJG(phi)
	IF (i<n-1) THEN
		t(i) = CONJG(ps(i))*ps(i)*gama(i+1)
		g = CONJG(phi)*gama(i+1)
 		END IF
	IF (i>1) THEN
		t(i-1) = CONJG(ps(i))*t(i-1)
	END IF
END DO
r(n) = aa
v(n) = one

!IF (MINVAL(ABS(r)) < eps*MAXVAL(ABS(r))) THEN
!	tr = 1.0d100	! tukaj je potrebno narediti stop, priblizek je dober ...
!	RETURN
!END IF

! Solve the linear system  Rw = v.
w(n) = v(n)/r(n)
w(n-1) = (v(n-1)-w(n)*s(n-1))/r(n-1)
DO i=n-2,1,-1
	w(i) = (v(i)-w(i+1)*s(i)-w(i+2)*t(i))/r(i)
END DO
	
! compute the trace
tr = SUM(Bf*u*w)
tr = tr - SUM(Bh*CONJG(ps)*u(1:n-1)*w(2:n))
tr = tr - SUM(Bg/CONJG(ps)*( u(2:n)*w(1:n-1) - (v(1:n-1)*u(2:n)+CONJG(ps(1:n-1))*ps(1:n-1))/r(1:n-1)))

END SUBROUTINE newtcorrQRcx

!======================================================================
!	newtcorrLUcx
!======================================================================

SUBROUTINE newtcorrLUcx(a, au, al, b, bu, bl, c, cu, cl, n, xval, tr)

! computes derivative of the determinant of a symmetric tridiagonal 
! quadratic matrix polynomial using LU and MU+LV decompositions
!	 
! f(x) = det(x^2*tridiag(al,a,au)+x*tridiag(bl,b,bu)+tridiag(cl,c,cu))
!
! --------- output    ------
! tr : f'(x)/f(x)
!
! --------- input     ------
! a,b,c : diagonal of A, B, C, respectively
! al, bl, cl : lower codiagonal of A, B, C, respectively
! au, bu, cu : upper codiagonal of A, B, C, respectively
! n : matrix dimension
! xval : real value
!
! --------- arguments ------
INTEGER						:: n
REAL(dp), DIMENSION(n)		:: a, b, c
REAL(dp), DIMENSION(n-1)	:: al, bl, cl, au, bu, cu
COMPLEX(dp)					:: xval, tr
!
! --------- local     ------
COMPLEX(dp), DIMENSION(n)		:: s, ds, u, x
COMPLEX(dp), DIMENSION(n-1)		:: l, tu, tl, dtu, dtl, v, y, m
COMPLEX(dp), DIMENSION(n-2)		:: w, z
COMPLEX(dp)						:: tmpu, tmpv
INTEGER							:: i, j, k, tmpp
INTEGER, DIMENSION(n)			:: p, inda, indb
    
s = xval**2*a + xval*b + c
tu = xval**2*au + xval*bu + cu
tl = xval**2*al + xval*bl + cl
ds = 2*xval*a + b
dtu = 2*xval*au + bu
dtl = 2*xval*al + bl

! LU decomposition with partial pivoting
! for the symmetric tridiagonal matrix tridiag(s,ds,ds)
	
DO i=1,n
	p(i) = i
END DO
inda = p-1
indb = inda
u(1) = s(1)
v(1) = tu(1)

z = zero
w = zero

DO k = 1,n-1
	! Pivoting 
	IF (ABS(u(k))<ABS(tl(k))) THEN
		tmpp = p(k+1)
		p(k+1) = p(k)
		p(k) = tmpp
        IF (k>1) THEN
			inda(k+1) = inda(k)
		END IF
		inda(k) = 0
		indb(k) = 0
		tmpu = u(k)
		tmpv = v(k)
		u(k) = tl(k)
		v(k) = s(k+1)
        l(k) = tmpu/tl(k)
        u(k+1) = tmpv-l(k)*v(k)
        IF (k<n-1) THEN
			w(k) = tu(k+1)
			v(k+1) = -l(k)*w(k)
        END IF
	ELSE
        l(k) = tl(k)/u(k)
        u(k+1) = s(k+1)-l(k)*v(k)
        IF (k<n-1) THEN
			v(k+1) = tu(k+1)
		END IF
    END IF
 END Do                

DO i = 1,n              
	! first part: elements of M are computed
    IF (inda(i) > 0) THEN
        DO j = inda(i),indb(i)
			IF (j == p(i)-1) THEN
				tmpu = dtl(j)    
			ELSE IF (j == p(i)) THEN
				tmpu = ds(j)
			ELSE IF (j == p(i) + 1) THEN
				tmpu = dtu(j-1)  
			ELSE
				tmpu = 0
			END IF	

			IF (inda(i) < j-1) THEN
				tmpu = tmpu - m(j-2)*w(j-2) - l(j-2)*z(j-2) - m(j-1)*v(j-1) - l(j-1)*y(j-1) - l(j)*x(j)
			ELSE IF (inda(i) == j) THEN
				tmpu = tmpu - l(j)*x(j)   
			ELSE IF (inda(i) == j-1) THEN
				tmpu = tmpu - m(j-1)*v(j-1) - l(j-1)*y(j-1) - l(j)*x(j)
			END IF
               
			m(j) = tmpu/u(j)
		END DO
	END IF
		
	! second part: elements of V (x,y,z)
	IF (i == p(i)-1) THEN
		tmpu = dtl(i)   
	ELSE IF (i == p(i)) THEN
		tmpu = ds(i)  
	ELSE IF (i == p(i)+1) THEN
		tmpu = dtu(i-1)  
	ELSE
		tmpu = 0
	END IF

	IF (inda(i)>0) THEN
		IF (inda(i)<i-1) THEN                       
			tmpu = tmpu - m(i-2)*w(i-2) - l(i-2)*z(i-2) - m(i-1)*v(i-1) - l(i-1)*y(i-1)
		ELSE
			tmpu = tmpu - m(i-1)*v(i-1) - l(i-1)*y(i-1)
		ENDIF
	END IF
        
	x(i) = tmpu
		
	IF (i<n) THEN
		IF (i+1 == p(i)) THEN
			tmpu = ds(i+1)
		ELSE IF (i+1 == p(i)+1) THEN
			tmpu = dtu(i)    
		ELSE IF (i+1 == p(i)-1) THEN
			tmpu = dtl(i+1)   
		ELSE
			tmpu = 0
		END IF

		IF (inda(i)>0) THEN
            tmpu = tmpu - m(i-1)*w(i-1) - l(i-1)*z(i-1)
		END IF	 
		
		y(i) = tmpu
	END IF

	IF (i<n-1) THEN
        IF (i+2==p(i)+1) THEN
			z(i) = dtu(i+1)   !!
		ELSE IF (i+2==p(i)) THEN
			z(i) = ds(i+2)
		ELSE
			tmpu = 0
		END IF
	END IF
END DO

tr = sum(x/u)

END SUBROUTINE newtcorrLUcx

!======================================================================
!	reigencx
!======================================================================

RECURSIVE SUBROUTINE reigencx(a, au, al, b, bu, bl, c, cu, cl, n, z, maxit, & 
						&   iterall, iterlast, itermax, mode)

! Recursive subroutine for the evaluation of the eigenvalues of a 
! tridiagonal hyperbolic quadratic eigenvalue problem 
! Q(x) = x^2*tridiag(al,a,au)+x*tridiag(bl,b,bu)+tridiag(cl,c,cu)
!
! --------- output    ------
! z : eigenvalues
! iterall  : all iterations
! iterlast : iterations in the last conquer phase
! itermax  : maximum number of iterations in recursive calls
!
! --------- input     ------
! a,b,c      : diagonal of A, B, C, respectively
! al, bl, cl : lower codiagonal of A, B, C, respectively
! au, bu, cu : upper codiagonal of A, B, C, respectively
! n          : matrix dimension
! maxit : maximum possible number of iterations
! mode  : method for the Newton correction: 1 (3-term), 2 (QR), 3 (LU)	
!
! --------- arguments ------
INTEGER                    :: n, maxit, iterall, iterlast, itermax, mode
REAL(dp), DIMENSION(n)     :: a, b, c
REAL(dp), DIMENSION(n-1)   :: al, bl, cl, au, bu, cu
COMPLEX(dp), DIMENSION(2*n):: z
!
! --------- local     ------
INTEGER                    :: i, nm, itermax1, itermax2, itermax3
INTEGER                    :: iter1, iter2, tmpiter
COMPLEX(dp)                :: discr, sqrtdisc
REAL                       :: rnd1, rnd2
REAL(dp)				   :: pertb	

! for n=1 the eigenvalues are solution of a quadratic equation
    
IF (n==1) THEN
    sqrtdisc = SQRT(CMPLX(b(1)**2-4*a(1)*c(1),0))
    IF (ABS(-b(1)-sqrtdisc)>=ABS(-b(1)+sqrtdisc)) THEN
		sqrtdisc = -sqrtdisc
	END IF
	z(1) = (- b(1) + sqrtdisc )/(2*a(1))
	Z(2) = c(1)/a(1)/z(1)
	iterall = 0
	iterlast = 0
	itermax = 0
    RETURN
 END IF
       
! problem of size n>1 is split into two subproblems

nm=(n+1)/2 
	   	
CALL reigencx(a(1:nm),au(1:nm-1),al(1:nm-1),b(1:nm),bu(1:nm-1),bl(1:nm-1),c(1:nm),cu(1:nm-1),cl(1:nm-1),& 
					& nm,z(1:2*nm),maxit,iter1,tmpiter,itermax1,mode)

CALL reigencx(a(nm+1:n),au(nm+1:n-1),al(nm+1:n-1),b(nm+1:n),bu(nm+1:n-1),bl(nm+1:n-1),c(nm+1:n),cu(nm+1:n-1),cl(nm+1:n-1),&
					& n-nm,z(2*nm+1:2*n),maxit,iter2,tmpiter,itermax2,mode)

! slightly perturb the eigenvalues of the two subproblems
! for avoiding multiple roots and merge them together

pertb = 10000*eps;
DO i=1,nm
	CALL RANDOM_NUMBER(rnd1)
	CALL RANDOM_NUMBER(rnd2)
	z(i) = z(i)*(one+pertb*CMPLX(0.5-rnd1,rnd2))
END DO
DO i=nm+1,n
	CALL RANDOM_NUMBER(rnd1)
	CALL RANDOM_NUMBER(rnd2)
	z(i) = z(i)*(one-pertb*CMPLX(0.5-rnd1,rnd2))
END DO

! refine the approximations by the Ehrlich-Aberth iteration
CALL aberthcx(a, au, al, b, bu, bl, c, cu, cl, n, z, maxit, iterlast, itermax3, mode)
iterall = iter1 + iter2 + iterlast   
itermax = MAX(itermax1,itermax2,itermax3)

END SUBROUTINE reigencx

!======================================================================
!	aberthcx
!======================================================================

SUBROUTINE aberthcx(a, au, al, b, bu, bl, c, cu, cl,n,z,maxit,iterall,iterlevel,mode)

! This subroutine performs the Aberth iterations for refining the 
! approximation z to the eigenvalue of
! Q(x) = x^2*tridiag(al,a,au)+x*tridiag(bl,b,bu)+tridiag(cl,c,cu)
!
! We use the following stopping criteria:
!	1) both the Newton and the Ehrlich-Aberth corrections are small
!	2) only some zeroes are left and corrections are not getting smaller
!
! --------- output    ------
! z : eigenvalues
! iterall : number of iterations required 
! iterlevel : maximum level of iteration
!
! --------- input     ------
! a,b,c      : diagonal of A, B, C, respectively
! al, bl, cl : lower codiagonal of A, B, C, respectively
! au, bu, cu : upper codiagonal of A, B, C, respectively
! n          : matrix dimension
! z			 : initial approximations
! maxit : maximum possible number of iterations
! mode  : method for the Newton correction: 1 (3-term), 2 (QR), 3 (LU)	
!
! --------- arguments ------
INTEGER                     :: n, maxit, iterall, iterlevel, mode
REAL(dp), DIMENSION(n)      :: a, b, c
REAL(dp), DIMENSION(n-1)	:: al, bl, cl, au, bu, cu
COMPLEX(dp), DIMENSION(2*n) :: z
!
! --------- local     ------
INTEGER                     :: i, j, m, nfound, razlog, lastchange
COMPLEX(dp)                 :: corr, abcorr, newt, den 
LOGICAL,DIMENSION(2*n)      :: cont
REAL(dp), DIMENSION(2*n)    :: minres
INTEGER, DIMENSION(2*n)		:: kand
REAL(dp)					:: relchange, rnd1, rnd2
INTEGER						:: mrem

cont = .TRUE.
nfound = 0
iterlevel = 0 
iterall = 0
m = 2*n
mrem = MIN((5*m)/10, m-2) ! number of eigenvalues for the 2nd criteria
kand = 0 ! number of consecutive iterations without residual improvement
minres = 1.d100 ! minimum residuals for the approximations

DO WHILE (nfound<m .AND. iterlevel<=maxit)
	iterlevel = iterlevel+1
    DO i=1,m
		IF (cont(i)) THEN
			iterall = iterall+1
            
			! Sum of reciprocal differences
			abcorr = 0
            DO j=1,m
				IF (.NOT. j == i) THEN
					IF (ABS(z(i)-z(j))==0) THEN
!						WRITE(11,*) "WARNING: zi=zj i=",i,"j=",j,"m=",m,"z(i)=",z(i)
						CALL RANDOM_NUMBER(rnd1)
						CALL RANDOM_NUMBER(rnd2)
						z(i) = z(i)*(one+100*eps*CMPLX(0.5-rnd1,rnd2))
					END IF
					abcorr = abcorr + one/(z(i)-z(j))
				END IF
            END DO

			! compute the Newton and the Aberth correction
			CALL newtcorrcx(a, au, al, b, bu, bl, c, cu, cl, n, z(i), newt, mode)
 			newt = one/newt
	            den = one-newt*abcorr
		        IF (den==0) THEN
			        corr = newt
				    WRITE(11,*)"WARNING: INF -> SWITCH TO Newton CORRECTION for i=",i, "n=",n, "z=",z(i), " corr=",corr
				ELSE
                    corr = newt/den
                END IF

		        IF (ISNAN(REAL(corr))) THEN
			        corr = 0
				    WRITE(11,*) "WARNING: NAN -> SWITCH TO ZERO CORRECTION for i=",i, "n=",n, "z=",z(i)
				END IF
        		z(i) = z(i)-corr
				relchange = abs(corr)/abs(z(i));
				IF (relchange>=minres(i)) THEN
					kand(i) = kand(i) + 1
				END IF
				razlog=0

			IF (cont(i)) THEN
				! First criteria : small correction
				IF ((ABS(newt)<epsrelaxed) .AND. (relchange<eps)) THEN
					razlog = 1
!				IF (relchange<eps) THEN
!					razlog = 1
				! Second criteria : no convergence for at least 5 cycles,
				! no residual improvement and the majority of the 
				! eigenvalues has already converged
				ELSE IF ((relchange<minres(i)*100) .AND. (kand(i)>2) &
				        & .AND. ((iterall-lastchange)>5*(m-nfound)) & 
				        & .AND. (nfound>=mrem) .AND. relchange<1.d-7) THEN
					razlog = 2
   				END IF
				IF (relchange<minres(i)) THEN
					minres(i) = relchange
	                kand(i) = 0
		        END IF
				IF (.NOT. razlog==0) THEN
					cont(i) = .FALSE.
					nfound = nfound + 1
					lastchange = iterall
				END IF
			END IF
		END IF  	   	
	END DO
END DO
IF (iterlevel>=maxit) THEN
   WRITE(11,*)"WARNING: Exceeded max iterations n=",n," iterlevel=",iterlevel," nfound=",nfound
!   STOP
END IF

END SUBROUTINE aberthcx

!======================================================================
!	newtcorrcx
!======================================================================

SUBROUTINE newtcorrcx(a, au, al, b, bu, bl, c, cu, cl,n,z,df,mode)

! Computes the Newton correction for the determinant of a real 
! symmetric tridiagonal quadratic polynomial matrix using
! 3-term, QR or LU
!
! --------- output    ------
! df : newton correction 
!
! --------- input     ------
! a,b,c      : diagonal of A, B, C, respectively
! al, bl, cl : lower codiagonal of A, B, C, respectively
! au, bu, cu : upper codiagonal of A, B, C, respectively
! n          : matrix dimension
! z			 : initial approximations
! mode  : method for the Newton correction: 1 (3-term), 2 (QR), 3 (LU)	
!
! --------- arguments ------
INTEGER						:: n, mode
REAL(dp), DIMENSION(n)      :: a, b, c
REAL(dp), DIMENSION(n-1)	:: al, bl, cl, au, bu, cu
COMPLEX(dp)					:: z, df

IF (mode==1) THEN
	CALL newtcorr3Tcx(a, au, al, b, bu, bl, c, cu, cl,n,z,df)  	! Three-term reccurences
ELSE IF (mode==2) THEN
	CALL newtcorrQRcx(a, au, al, b, bu, bl, c, cu, cl,n,z,df)	! QR approach
ELSE IF (mode==3) THEN
	CALL newtcorrLUcx(a, au, al, b, bu, bl, c, cu, cl,n,z,df)	! LU (MU+LV) approach
ELSE
	df = zero
END IF

END SUBROUTINE newtcorrcx

!======================================================================

END MODULE qep3deacx

