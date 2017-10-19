PROGRAM conv_test
IMPLICIT NONE
INTEGER                                         :: deg, itmax, test_num
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: p, berr, er, ei, error
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: exacteigs
!loop variables
INTEGER                                         :: i, j
!intrinsic subroutines
INTRINSIC                                       :: dble, dcmplx, dimag, epsilon
!external subroutines
EXTERNAL                                        :: dslm_conv, dsam_conv, dsam
!array varaiables
CHARACTER(LEN=100)                              :: arg

itmax = 60
test_num = 127 

CALL getarg(1,arg)  
READ(arg, *) deg

 OPEN(UNIT=1,FILE="results/dslm_conv_test_results.txt")
DO j=1,test_num
  !Create random polynomial
  ALLOCATE(p(deg+1))
  CALL daruv(deg+1,p)

  !DSLM 
  ALLOCATE(berr(deg),er(deg),ei(deg), exacteigs(deg), error(itmax+1))

  !Compute convergent vector
  CALL dslm(p, deg, er, ei, berr)
  DO i=1,deg
    exacteigs(i) = dcmplx(er(i), ei(i))
  ENDDO
  DEALLOCATE(berr,er,ei)

  !Test for convergence rate
  ALLOCATE(berr(deg),er(deg),ei(deg))
  CALL dslm_conv(p, deg, er, ei, berr, exacteigs, error)
  DO i =1,itmax
    WRITE(1, '(ES15.2E3)') error(i)  
  ENDDO
  DEALLOCATE(berr,er,ei, exacteigs, error, p)

  WRITE(1,*)

ENDDO

 CLOSE(UNIT=1)

CALL EXECUTE_COMMAND_LINE('python conv_test.py')

END PROGRAM conv_test
