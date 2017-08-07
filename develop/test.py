from ctypes import *
from random import *
import sys
lib = cdll.LoadLibrary(“./libtest.a”)

startDegree=sys.argv[1]
maxDegree=sys.argv[2]


degree=startDegree
while degree<maxDegree:
	

	method = libadd.dseval
	#
	# The byref() is necessary since
	# FORTRAN does references,
	# and not values (like e.g. C)
	#
	method( byref(x), byref(y) )


