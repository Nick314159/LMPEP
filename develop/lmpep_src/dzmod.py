"""
author Thomas R. Cameron, Davidson College
author Nikolas I. Steckley, Portland State University
date 2017
brief <b> DZMOD returns moduli of complex number. </b>
par Purpose:
verbatim
!> DZMOD returns moduli of complex number a+bi, while avoiding harmul overflow and underflow. 
endverbatim
param[in] a
verbatim Double precision number, real part. \endverbatim
param[in] b
verbatim Double precision number, imaginary part.\endverbatim
"""
import math
epsilon = sys.float_info.epsilon
zero = 0.0

def dzmod(a, b):
    if abs(a) < epsilon and abs(b) < epsilon:
        dzmod = zero
    elif abs(a) < abs(b):
        dzmod = abs(b) * sqrt(1+(a/b)**2)
    else:
        dzmod = abs(a) * sqrt(1+(b/a)**2)
    return dzmod


