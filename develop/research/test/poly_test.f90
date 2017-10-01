PROGRAM poly_test
IMPLICIT NONE
INTEGER                                         :: deg
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: p
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: exacteigs
DOUBLE PRECISION, PARAMETER                     :: pi = 3.1415926535897932d0
!intrinsic subroutines
INTRINSIC                                       :: cos, dcmplx

OPEN(UNIT=1,FILE="results/poly_test_results.txt")

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

!Rev Wilkinson Deg 10 Polynomial
deg=10
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=1D0/3628800D0
p(2)=-11D0/725760D0
p(3)=11D0/30240D0
p(4)=-121D0/24192D0
p(5)=7513D0/172800D0
p(6)=-8591D0/34560D0
p(7)=341693D0/362880D0
p(8)=-84095D0/36288D0
p(9)=177133D0/50400D0
p(10)=-7381D0/2520D0
p(11)=1D0
exacteigs(1)=dcmplx(1D0,0D0)
exacteigs(2)=dcmplx(1D0/2D0,0D0)
exacteigs(3)=dcmplx(1D0/3D0,0D0)
exacteigs(4)=dcmplx(1D0/4D0,0D0)
exacteigs(5)=dcmplx(1D0/5D0,0D0)
exacteigs(6)=dcmplx(1D0/6D0,0D0)
exacteigs(7)=dcmplx(1D0/7D0,0D0)
exacteigs(8)=dcmplx(1D0/8D0,0D0)
exacteigs(9)=dcmplx(1D0/9D0,0D0)
exacteigs(10)=dcmplx(1D0/10D0,0D0)
WRITE(1, '(A)')  "Deg 10 Rev Wilkins Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

!Rev Wilkinson Deg 15 Polynomial
deg=15
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=-1D0/1307674368000D0
p(2)=1D0/10897286400D0
p(3)=-47D0/9340531200D0
p(4)=1D0/5987520D0
p(5)=-26921D0/7185024000D0
p(6)=109D0/1814400D0
p(7)=-324509D0/457228800D0
p(8)=4783D0/762048D0
p(9)=-54576553D0/1306368000D0
p(10)=2271089D0/10886400D0
p(11)=-277382447D0/359251200D0
p(12)=2065639D0/997920D0
p(13)=-35118025721D0/9081072000D0
p(14)=13215487D0/2802800D0
p(15)=-1195757D0/360360D0
p(16)=1D0
exacteigs(1)=dcmplx(1D0,0D0)
exacteigs(2)=dcmplx(1D0/2D0,0D0)
exacteigs(3)=dcmplx(1D0/3D0,0D0)
exacteigs(4)=dcmplx(1D0/4D0,0D0)
exacteigs(5)=dcmplx(1D0/5D0,0D0)
exacteigs(6)=dcmplx(1D0/6D0,0D0)
exacteigs(7)=dcmplx(1D0/7D0,0D0)
exacteigs(8)=dcmplx(1D0/8D0,0D0)
exacteigs(9)=dcmplx(1D0/9D0,0D0)
exacteigs(10)=dcmplx(1D0/10D0,0D0)
exacteigs(11)=dcmplx(1D0/11D0,0D0)
exacteigs(12)=dcmplx(1D0/12D0,0D0)
exacteigs(13)=dcmplx(1D0/13D0,0D0)
exacteigs(14)=dcmplx(1D0/14D0,0D0)
exacteigs(15)=dcmplx(1D0/15D0,0D0)
WRITE(1, '(A)') "Deg 15 Rev Wilkins Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

!Rev Wilkinson Deg 20 Polynomial
!deg=20
!ALLOCATE(exacteigs(deg), p(deg+1))
!p(1)=
!p(2)=
!p(3)=
!p(4)=
!p(5)=
!p(6)=
!p(7)=
!p(8)=
!p(9)=
!p(10)=
!p(11)=
!p(12)=
!p(13)=
!p(14)=
!p(15)=
!p(16)=
!p(17)=
!p(18)=
!p(19)=
!p(20)=
!p(21)=1D0

!exacteigs(1)=dcmplx(1D0,0D0)
!exacteigs(2)=dcmplx(1D0/2D0,0D0)
!exacteigs(3)=dcmplx(1D0/3D0,0D0)
!exacteigs(4)=dcmplx(1D0/4D0,0D0)
!exacteigs(5)=dcmplx(1D0/5D0,0D0)
!exacteigs(6)=dcmplx(1D0/6D0,0D0)
!exacteigs(7)=dcmplx(1D0/7D0,0D0)
!exacteigs(8)=dcmplx(1D0/8D0,0D0)
!exacteigs(9)=dcmplx(1D0/9D0,0D0)
!exacteigs(10)=dcmplx(1D0/10D0,0D0)
!exacteigs(11)=dcmplx(1D0/11D0,0D0)
!exacteigs(12)=dcmplx(1D0/12D0,0D0)
!exacteigs(13)=dcmplx(1D0/13D0,0D0)
!exacteigs(14)=dcmplx(1D0/14D0,0D0)
!exacteigs(15)=dcmplx(1D0/15D0,0D0)
!exacteigs(16)=dcmplx(1D0/16D0,0D0)
!exacteigs(17)=dcmplx(1D0/17D0,0D0)
!exacteigs(18)=dcmplx(1D0/18D0,0D0)
!exacteigs(19)=dcmplx(1D0/19D0,0D0)
!exacteigs(20)=dcmplx(1D0/20D0,0D0)
!WRITE(1, '(A)')  "Deg 20 Rev Wilkins Poly"
!CALL test(deg, p, exacteigs)
!DEALLOCATE(exacteigs, p)

 CLOSE(UNIT=1)

CALL EXECUTE_COMMAND_LINE('python poly_test.py')

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
SUBROUTINE zsort(root, n)
IMPLICIT NONE
scalar arguments
INTEGER, INTENT(IN)             :: n
array arguments
DOUBLE COMPLEX, INTENT(INOUT)   :: root(*)
local scalars
INTEGER                         :: i, j, k
DOUBLE COMPLEX                  :: temp

DO i=1,n
    temp=root(i); j=i
    DO k=i+1,n
        IF (dble(temp)>dble(root(k))) THEN
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
