import csv
import math
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
npRealEstimate, npImaginaryEstimate, nrRealEstimate, nrImaginaryEstimate, nrReal, nrImaginary,  = [], [], [], [], [], [], 
#npReal, npImaginary = [], []
tests=['butterfly','cd_player','spring']

def format(s):
  str = s.strip()
  f = float(str)
  r = round(f, 6)
  return r

for test in tests :
  with open('../results/outputIepoly2-'+test+'.txt') as f:
    reader = csv.reader(f)
    reader = list(reader)
    row_count = sum(1 for row in reader)
    row_index = 2
    while True:
      row = reader[row_index]
      row_index = row_index + 1
      if not row : break
      nrRealEstimate.append(format(row[0]))
      nrImaginaryEstimate.append(format(row[1]))	
    while True:
      row = reader[row_index]
      row_index = row_index + 1
      if not row : break
      nrReal.append(format(row[0]))
      nrImaginary.append(format(row[1]))
    
    row_index = row_index + 2
    while True:
      row = reader[row_index]
      row_index = row_index + 1
      if not row : break
      npRealEstimate.append(format(row[0]))
      npImaginaryEstimate.append(format(row[1]))
 
  fig = plt.figure()
  ax = fig.add_subplot(111)
  
  ax.plot(npRealEstimate, npImaginaryEstimate, 'wo', label='NP Initial Estimate')
  ax.plot(nrRealEstimate, nrImaginaryEstimate, 'ko', label='NR Initial Estimate')
  ax.plot(nrReal, nrImaginary, 'k*', label='Eigenvalues')
  ax.set_ylabel('Real')
  ax.set_xlabel('Imaginary')
  legend = ax.legend(loc=0, shadow=True)
  savefig("../results/iepoly2"+test+".pdf")


