import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
npRealEstimate, npImaginaryEstimate, nrRealEstimate, nrImaginaryEstimate, npReal, npImaginary, nrReal, nrImaginary = [], [], [], [], [], [], [], []
tests=['bilby','butterfly','cd_player','spring']

for test in tests :
  with open('../results/outputIepoly2-'+test+'.txt') as f:
    print(test)
    reader = csv.reader(f)
    reader = list(reader)
    row_count = sum(1 for row in reader)
    row_index = 2
    while True:
      row = reader[row_index]
      row_index = row_index + 1
      if not row : break
      nrRealEstimate.append(float(row[0].strip()))
      nrImaginaryEstimate.append(float(row[1].strip()))	
    while True:
      row = reader[row_index]
      row_index = row_index + 1
      if not row : break
      nrReal.append(float(row[0].strip()))
      nrImaginary.append(float(row[1].strip()))
    
    row_index = row_index + 2
    while True:
      row = reader[row_index]
      row_index = row_index + 1
      if not row : break
      npRealEstimate.append(float(row[0].strip()))
      npImaginaryEstimate.append(float(row[1].strip()))
    while True:
      row = reader[row_index]
      row_index = row_index + 1
      if not row : break
      npReal.append(float(row[0].strip()))
      npImaginary.append(float(row[1].strip()))
      
      
  fig, ax = plt.subplots()
  
  maxReal = max(max(nrReal), max(nrRealEstimate), max(npReal), max(npRealEstimate))
  maxImaginary = max(max(nrImaginary), max(nrImaginaryEstimate), max(npImaginary), max(npImaginaryEstimate))
  minReal = min(min(nrReal), min(nrRealEstimate), min(npReal), min(npRealEstimate))
  minImaginary = min(max(nrImaginary), min(nrImaginaryEstimate), min(npImaginary), min(npImaginaryEstimate))
  ax.set_xlim(minReal, maxReal)
  ax.set_ylim(minImaginary, maxImaginary)
  ax.scatter(nrRealEstimate, nrImaginaryEstimate, c='r', marker='o', label='NR Initial Estimate')
  ax.scatter(nrReal, nrImaginary, c='r', marker='.', label='NR Actual')
  ax.scatter(npRealEstimate, npImaginaryEstimate, c='b', marker='o', label='NP Initial Estimate')
  ax.scatter(npReal, npImaginary, c='b', marker='.', label='NP Actual')
  ax.set_ylabel('Real')
  ax.set_xlabel('Imaginary')
  legend = ax.legend(loc=0, shadow=True)
  savefig("../results/iepoly2"+test+".pdf")


