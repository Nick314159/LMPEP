from ctypes import *
from random import *
import time
import sys
#Import Fortran library
lib = cdll.LoadLibrary("/libtest.a")

#Make methods available for execution
dseval = lib.dseval
drarr = lib.drarr
drnum = lib.drnum

zseval = lib.zseval
zrarr = lib.zrarr
zrnum = lib.zrnum

#Command line arguements
startDegree=sys.argv[1]
maxDegree=sys.argv[2]
dt=sys.argv[3]

#Hard coded values
itmax =100
der = 2

t, a = 0 #Need to allocate memory for this

degree=startDegree
if dt=="D":
	while degree<maxDegree:
		p = [None]* (degree+1)
		start = time.clock
		for i in range(itmax):
			# NOTE: The byref() is necessary since
			# FORTRAN does references,
			# and not values (like e.g. C)
			drarr(byref(p), byref(degree)+1)
			drnum(byref(t))
			dseval(byref(p), byref(t), byref(degree), byref(der), byref(a))
		endTime = clock.time
		print("Degree, Average Time")
		print(degree + "," + (startTime-endTime)/itmax)
		degree *= 2
	

elif dt=="Z":
	print(1)
