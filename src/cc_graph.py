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
    degree.append(float(row[0].strip()))
    getimeDegree.append(float(row[2].strip()))
    hstimeDegree.append(float(row[3].strip()))
with open('../results/outputComplexitySize.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    size.append(float(row[1].strip()))
    getimeSize.append(float(row[2].strip()))
    hstimeSize.append(float(row[3].strip()))

degree=np.array(degree)
getimeDegree=np.array(getimeDegree)

fig, ax = plt.subplots()
ax.loglog(degree, getimeDegree, 'k-o', label='General')
ax.loglog(degree, getimeDegree[0]*(degree/degree[0])**2, 'k--', label='Quadratic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_ge_times_degree.pdf")

fig, ax = plt.subplots()
ax.semilogy(degree, hstimeDegree, 'k-o', label='Hessenberg')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_hs_times_degree.pdf")

fig, ax = plt.subplots()
ax.semilogy(size, getimeSize, 'k-o', label='General')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_ge_times_size.pdf")

fig, ax = plt.subplots()
ax.semilogy(size, hstimeSize, 'k-o', label='Hessenberg')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_hs_times_size.pdf")


