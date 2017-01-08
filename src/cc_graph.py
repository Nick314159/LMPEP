import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, size, getimeDegree, hstimeDegree, getimeSize, hstimeSize = [], [], [], [], [], []

with open('../results/outputComplexityDegree.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(row[0].strip())
    getimeDegree.append(row[2].strip())
    hstimeDegree.append(row[3].strip())
with open('../results/outputComplexitySize.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    size.append(row[1].strip())
    getimeSize.append(row[2].strip())
    hstimeSize.append(row[3].strip())

fig, ax = plt.subplots()
ax.semilogy(degree, getimeDegree, 'k-*', label='General')
ax.semilogy(degree, hstimeDegree, 'k--', label='Hessenburg')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_times_degree.pdf")

fig, ax = plt.subplots()
ax.semilogy(size, getimeSize, 'k-*', label='General')
ax.semilogy(size, hstimeSize, 'k--', label='Hessenburg')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_times_size.pdf")

