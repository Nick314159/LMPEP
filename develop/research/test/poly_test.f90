PROGRAM poly_test
IMPLICIT NONE
INTEGER                                         :: deg
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: p
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: exacteigs
DOUBLE PRECISION, PARAMETER                     :: pi = 3.1415926535897932d0
!intrinsic subroutines
INTRINSIC                                       :: cos, dcmplx

OPEN(UNIT=1,FILE="poly_test_results.txt")

!Chebysehv Deg 20 Polynomial
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(21) = 1d0
p(20) = 0d0
p(19) =-2621440d0/524288d0
p(18) = 0d0
p(17) =+5570560d0/524288d0
p(16) = 0d0
p(15) =-6553600d0/524288d0
p(14) = 0d0
p(13) =+4659200d0/524288d0
p(12) = 0d0
p(11)=-2050048d0/524288d0
p(10)= 0d0
p(9)=+549120d0/524288d0
p(8)= 0d0
p(7)=-84480d0/524288d0
p(6)= 0d0
p(5)=+6600d0/524288d0
p(4)= 0d0
p(3)=-200d0/524288d0
p(2)= 0d0
p(1)= 1d0/524288d0
exacteigs(20)  = dcmplx(cos(1d0/40d0*pi),0d0)
exacteigs(19)  = dcmplx(cos(3d0/40d0*pi),0d0)
exacteigs(18)  = dcmplx(cos(5d0/40d0*pi),0d0)
exacteigs(17)  = dcmplx(cos(7d0/40d0*pi),0d0)
exacteigs(16)  = dcmplx(cos(9d0/40d0*pi),0d0)
exacteigs(15)  = dcmplx(cos(11d0/40d0*pi),0d0)
exacteigs(14)  = dcmplx(cos(13d0/40d0*pi),0d0)
exacteigs(13)  = dcmplx(cos(15d0/40d0*pi),0d0)
exacteigs(12)  = dcmplx(cos(17d0/40d0*pi),0d0)
exacteigs(11)  = dcmplx(cos(19d0/40d0*pi),0d0)
exacteigs(10)  = dcmplx(cos(21d0/40d0*pi),0d0)
exacteigs(9)  = dcmplx(cos(23d0/40d0*pi),0d0)
exacteigs(8)  = dcmplx(cos(25d0/40d0*pi),0d0)
exacteigs(7)  = dcmplx(cos(27d0/40d0*pi),0d0)
exacteigs(6)  = dcmplx(cos(29d0/40d0*pi),0d0)
exacteigs(5)  = dcmplx(cos(31d0/40d0*pi),0d0)
exacteigs(4)  = dcmplx(cos(33d0/40d0*pi),0d0)
exacteigs(3)  = dcmplx(cos(35d0/40d0*pi),0d0)
exacteigs(2)  = dcmplx(cos(37d0/40d0*pi),0d0)
exacteigs(1)  = dcmplx(cos(39d0/40d0*pi),0d0)
WRITE(1, '(A)') "Deg 20 Chebyshev Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

!Wilkinson Deg 10 Polynomial
deg=10
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=3628800D0
p(2)=-10628640D0
p(3)=12753576D0
p(4)=-8409500D0
p(5)=3416930D0
p(6)=-902055D0
p(7)=157773D0
p(8)=-18150D0
p(9)=1320D0
p(10)=-55D0
p(11)=1D0
exacteigs(1)=dcmplx(1D0,0D0)
exacteigs(2)=dcmplx(2D0,0D0)
exacteigs(3)=dcmplx(3D0,0D0)
exacteigs(4)=dcmplx(4D0,0D0)
exacteigs(5)=dcmplx(5D0,0D0)
exacteigs(6)=dcmplx(6D0,0D0)
exacteigs(7)=dcmplx(7D0,0D0)
exacteigs(8)=dcmplx(8D0,0D0)
exacteigs(9)=dcmplx(9D0,0D0)
exacteigs(10)=dcmplx(10D0,0D0)
WRITE(1, '(A)')  "Deg 10 Wilkins Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

!Wilkinson Deg 15 Polynomial
deg=15
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=-1307674368000D0
p(2)=4339163001600D0
p(3)=-6165817614720D0
p(4)=5056995703824D0
p(5)=-2706813345600D0
p(6)=1009672107080D0
p(7)=-272803210680D0
p(8)=54631129553D0
p(9)=-8207628000D0
p(10)=928095740D0
p(11)=-78558480D0
p(12)=4899622D0
p(13)=-218400D0
p(14)=6580D0
p(15)=-120D0
p(16)=1D0
exacteigs(1)=dcmplx(1D0,0D0)
exacteigs(2)=dcmplx(2D0,0D0)
exacteigs(3)=dcmplx(3D0,0D0)
exacteigs(4)=dcmplx(4D0,0D0)
exacteigs(5)=dcmplx(5D0,0D0)
exacteigs(6)=dcmplx(6D0,0D0)
exacteigs(7)=dcmplx(7D0,0D0)
exacteigs(8)=dcmplx(8D0,0D0)
exacteigs(9)=dcmplx(9D0,0D0)
exacteigs(10)=dcmplx(10D0,0D0)
exacteigs(11)=dcmplx(11D0,0D0)
exacteigs(12)=dcmplx(12D0,0D0)
exacteigs(13)=dcmplx(13D0,0D0)
exacteigs(14)=dcmplx(14D0,0D0)
exacteigs(15)=dcmplx(15D0,0D0)
WRITE(1, '(A)') "Deg 15 Wilkins Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

!Wilkinson Deg 20 Polynomial
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=2432902008176640000D0
p(2)=-8752948036761600000D0
p(3)=13803759753640704000D0
p(4)=-12870931245150988800D0
p(5)=8037811822645051776D0
p(6)=-3599979517947607200D0
p(7)=1206647803780373360D0
p(8)=-311333643161390640D0
p(9)=63030812099294896D0
p(10)=-10142299865511450D0
p(11)=1307535010540395D0
p(12)=-135585182899530D0
p(13)=11310276995381D0
p(14)=-756111184500D0
p(15)=40171771630D0
p(16)=-1672280820D0
p(17)=53327946D0
p(18)=-1256850D0
p(19)=20615D0
p(20)=-210D0
p(21)=1D0

exacteigs(1)=dcmplx(1D0,0D0)
exacteigs(2)=dcmplx(2D0,0D0)
exacteigs(3)=dcmplx(3D0,0D0)
exacteigs(4)=dcmplx(4D0,0D0)
exacteigs(5)=dcmplx(5D0,0D0)
exacteigs(6)=dcmplx(6D0,0D0)
exacteigs(7)=dcmplx(7D0,0D0)
exacteigs(8)=dcmplx(8D0,0D0)
exacteigs(9)=dcmplx(9D0,0D0)
exacteigs(10)=dcmplx(10D0,0D0)
exacteigs(11)=dcmplx(11D0,0D0)
exacteigs(12)=dcmplx(12D0,0D0)
exacteigs(13)=dcmplx(13D0,0D0)
exacteigs(14)=dcmplx(14D0,0D0)
exacteigs(15)=dcmplx(15D0,0D0)
exacteigs(16)=dcmplx(16D0,0D0)
exacteigs(17)=dcmplx(17D0,0D0)
exacteigs(18)=dcmplx(18D0,0D0)
exacteigs(19)=dcmplx(19D0,0D0)
exacteigs(20)=dcmplx(20D0,0D0)
WRITE(1, '(A)')  "Deg 20 Wilkins Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

!Cylcotomic Deg 20 Polynomial
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=1D0
p(2)=0D0
p(3)=0D0
p(4)=0D0
p(5)=0D0
p(6)=1D0
p(7)=0D0
p(8)=0D0
p(9)=0D0
p(10)=0D0
p(11)=1D0
p(12)=0D0
p(13)=0D0
p(14)=0D0
p(15)=0D0
p(16)=1D0
p(17)=0D0
p(18)=0D0
p(19)=0D0
p(20)=0D0
p(21)=1D0
exacteigs(1)=dcmplx(-0.992114701314478, -0.125333233564304)
exacteigs(2)=dcmplx(0.968583161128631, 0.248689887164855)
exacteigs(3)=dcmplx(-0.929776485888251, -0.368124552684678)
exacteigs(4)=dcmplx(0.876306680043864, 0.481753674101715)
exacteigs(5)=dcmplx(0.728968627421412, 0.684547105928689)
exacteigs(6)=dcmplx(-0.637423989748690, -0.770513242775789)
exacteigs(7)=dcmplx(0.535826794978997, 0.844327925502015)
exacteigs(8)=dcmplx(-0.425779291565073, -0.904827052466020)
exacteigs(9)=dcmplx(-0.187381314585725, -0.982287250728689)
exacteigs(10)=dcmplx(0.062790519529313, 0.998026728428272)
exacteigs(11)=dcmplx(0.062790519529313, -0.998026728428272)
exacteigs(12)=dcmplx(-0.187381314585725, 0.982287250728689)
exacteigs(13)=dcmplx(-0.425779291565073, 0.904827052466020)
exacteigs(14)=dcmplx(0.535826794978997, -0.844327925502015)
exacteigs(15)=dcmplx(-0.637423989748690, 0.770513242775789)
exacteigs(16)=dcmplx(0.728968627421412, -0.684547105928689)
exacteigs(17)=dcmplx(0.876306680043864, -0.481753674101715)
exacteigs(18)=dcmplx(-0.929776485888251, 0.368124552684678)
exacteigs(19)=dcmplx(0.968583161128631, -0.248689887164855)
exacteigs(20)=dcmplx(-0.992114701314478, 0.125333233564304)
WRITE(1, '(A)')  "Deg 20 Cyclotomic Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

 CLOSE(UNIT=1)
CONTAINS

!********************************************************
!                       TEST                            *
!******************************************************** 
SUBROUTINE test(deg, p, exacteigs)
USE poly_zeroes
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                             :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)                    :: p(*) 
DOUBLE COMPLEX, INTENT(IN)                      :: exacteigs(*)
!LMPEP variables
DOUBLE PRECISION, DIMENSION(deg)                :: berr, er, ei
!POLZEROS variables
INTEGER                                         :: iter, nitmax
LOGICAL, DIMENSION(deg+1)                       :: err
REAL(KIND=dp), DIMENSION(1:deg)                 :: radius
COMPLEX(KIND=dp), DIMENSION(0:deg)              :: poly
COMPLEX(KIND=dp), DIMENSION(1:deg)              :: root
REAL(KIND=dp), PARAMETER                        :: eps=epsilon(1.0D0), big=huge(1.0D0), small=tiny(1.0D0)
!AMVW variables
INTEGER                                         :: FLAG
INTEGER, DIMENSION(deg)                         :: ITS
DOUBLE PRECISION, DIMENSION(deg)                :: REIGS, IEIGS, P2
!loop variables
INTEGER                                         :: i
!intrinsic subroutines
INTRINSIC                                       :: dble, dcmplx, dimag, epsilon, tiny, huge
!external subroutines
EXTERNAL                                        :: dslm, damvw
!external functions
DOUBLE PRECISION                                :: dzmod
EXTERNAL                                        :: dzmod


!LMPEP
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(A)') 'LMPEP Absolute Error:'
DO i=1,deg
   WRITE(1, *) dzmod(dble(exacteigs(i))-er(i),dimag(exacteigs(i))-ei(i))
ENDDO
PRINT*, ''

!POLZEROS
nitmax=60
DO i=0,deg
    poly(i)=dcmplx(p(i+1), 0.0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(A)') 'POLZEROS Absolute Error:'
DO i=1,deg
    WRITE(1, *) dzmod(dble(exacteigs(i))-dble(root(i)),dimag(exacteigs(i))-dimag(root(i)))
ENDDO
PRINT*, ''

!AMVW
FLAG=1
DO i=1,deg
    p2(deg-i+1)=p(i)/p(deg+1)
ENDDO
CALL damvw(deg, P2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS, IEIGS, deg)
WRITE(1, '(A)') 'AMVW Absolute Error:'
DO i=1,deg
    WRITE(1, *) dzmod(dble(exacteigs(i))-REIGS(i),dimag(exacteigs(i))-IEIGS(i))
ENDDO
PRINT*, ''

END SUBROUTINE test

!********************************************************
!                       DSORT                           *
!******************************************************** 
SUBROUTINE dsort(er, ei, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)             :: n
!array arguments
DOUBLE PRECISION, INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER                         :: i, j, k
DOUBLE PRECISION                :: tr, ti

DO i=1,n
    tr=er(i); ti=ei(i); j=i
    DO k=i+1,n
        IF (tr>er(k)) THEN
            tr=er(k); ti=ei(k); j=k
        ENDIF
    ENDDO
    IF (j>i) THEN
        tr=er(i); ti=ei(i)
        er(i)=er(j); ei(i)=ei(j)
        er(j)=tr; ei(j)=ti
    ENDIF
ENDDO
END SUBROUTINE dsort

!********************************************************
!                       ZSORT                           *
!********************************************************             
!SUBROUTINE zsort(root, n)
!IMPLICIT NONE
!scalar arguments
!INTEGER, INTENT(IN)             :: n
!array arguments
!DOUBLE COMPLEX, INTENT(INOUT)   :: root(*)
!local scalars
!INTEGER                         :: i, j, k
!DOUBLE COMPLEX                  :: temp

!DO i=1,n
!    temp=root(i); j=i
!    DO k=i+1,n
!        IF (dble(temp)>dble(root(k))) THEN
!            temp=root(k); j=k
!        ENDIF
!    ENDDO
!    IF (j>i) THEN
!        temp=root(i)
!        root(i)=root(j)
!        root(j)=temp
!    ENDIF
!ENDDO
!END SUBROUTINE zsort

!********************************************************
!                       ZSORT                           *
!********************************************************             
SUBROUTINE zsort(root, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)             :: n
!array arguments
DOUBLE COMPLEX, INTENT(INOUT)   :: root(*)
!local scalars
INTEGER                         :: i, j, k
DOUBLE COMPLEX                  :: temp

DO i=1,n
    temp=root(i); j=i
    DO k=i+1,n
        IF (atan2(aimag(temp), real(temp))>atan2(aimag(root(k)), real(root(k)))) THEN
            temp=root(k); j=k
        ENDIF
    ENDDO
    IF (j>i) THEN
        temp=root(i)
        root(i)=root(j)
        root(j)=temp
    ENDIF
ENDDO
END SUBROUTINE zsort
END PROGRAM poly_test