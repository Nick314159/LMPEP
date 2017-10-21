PROGRAM conv_test
IMPLICIT NONE
INTEGER                                         :: deg, i
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: p
DOUBLE PRECISION, PARAMETER                     :: pi = 3.1415926535897932d0

!intrinsic subroutines
INTRINSIC                                       :: cos, dcmplx

OPEN(UNIT=1,FILE="results/dslm_conv_test_results.txt")

!Chebysehv Deg 20 Polynomial
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)') "Deg 20 Chebyshev Poly"
CALL test(deg, p)
DEALLOCATE(p)

!Wilkinson Deg 10 Polynomial
deg=10
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)')  "Deg 10 Wilkins Poly"
CALL test(deg, p)
DEALLOCATE(p)


!Wilkinson Deg 15 Polynomial
deg=15
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)') "Deg 15 Wilkins Poly"
CALL test(deg, p)
DEALLOCATE(p)

!Wilkinson Deg 20 Polynomial
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)')  "Deg 20 Wilkins Poly"
CALL test(deg, p)
DEALLOCATE(p)

!Rev Wilkinson Deg 10 Polynomial
deg=10
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)')  "Deg 10 Rev Wilkins Poly"
CALL test(deg, p)
DEALLOCATE(p)

!Rev Wilkinson Deg 15 Polynomial
deg=15
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)') "Deg 15 Rev Wilkins Poly"
CALL test(deg, p)
DEALLOCATE(p)


!Create random polynomial
deg=20
ALLOCATE(p(deg+1))
CALL daruv(deg+1,p)
WRITE(1, '(A)') "Deg 20 Random Poly"
CALL test(deg, p)
DEALLOCATE(p)

 CLOSE(UNIT=1)

CALL EXECUTE_COMMAND_LINE('python conv_test.py')

CONTAINS

!********************************************************
!                       TEST                            *
!******************************************************** 
SUBROUTINE test(deg, p)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                             :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)                    :: p(*) 
!LMPEP variables
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: berr, er, ei, exacteigs, error
!loop variables
INTEGER                                         :: i
!intrinsic subroutines
INTRINSIC                                       :: dble, dimag, dcmplx
!external subroutines
EXTERNAL                                        :: dslm
!external functions
DOUBLE PRECISION                                :: dzmod
EXTERNAL                                        :: dzmod
INTEGER, PARAMETER                              :: itmax = 60

ALLOCATE(berr(deg),er(deg),ei(deg), exacteigs(deg), error(itmax+1))
CALL dslm(p, deg, er, ei, berr)
DO i=1,deg
   exacteigs=dcmplx(er(i), ei(i))
ENDDO
DEALLOCATE(er, ei, berr)

!Test for convergence rate
ALLOCATE(berr(deg),er(deg),ei(deg))
CALL dslm_conv(p, deg, er, ei, berr, exacteigs, error)
DO i=1,itmax+1
  WRITE(1, '(I10)', advance='no') i
  WRITE(1, '(A)', advance='no') ','
  WRITE(1, '(ES15.2)') error(i)
ENDDO
DEALLOCATE(berr,er,ei, exacteigs, error)

END SUBROUTINE test

END PROGRAM conv_test
