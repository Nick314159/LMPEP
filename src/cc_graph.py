import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, size, getimeDegree, hstimeDegree, tdtimeDegree, getimeSize, hstimeSize, tdtimeSize = [], [], [], [], [], [], [], []

with open('../results/outputComplexityDegree.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(float(row[0].strip()))
    getimeDegree.append(float(row[2].strip()))
    hstimeDegree.append(float(row[3].strip()))
    tdtimeDegree.append(float(row[4].strip()))
with open('../results/outputComplexitySize.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    size.append(float(row[1].strip()))
    getimeSize.append(float(row[2].strip()))
    hstimeSize.append(float(row[3].strip()))
    tdtimeSize.append(float(row[4].strip()))

degree=np.array(degree)
size=np.array(size)
getimeDegree=np.array(getimeDegree)
hstimeDegree=np.array(hstimeDegree)
tdtimeDegree=np.array(tdtimeDegree)
getimeSize=np.array(getimeSize)
hstimeSize=np.array(hstimeSize)
tdtimeSize=np.array(tdtimeSize)


fig, ax = plt.subplots()
ax.loglog(degree, getimeDegree, 'k-o', label='General')
ax.loglog(degree, getimeDegree[0]*(degree/degree[0])**2, 'k--', label='Quadratic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_ge_times_degree.pdf", format='pdf')

fig, ax = plt.subplots()
ax.loglog(degree, hstimeDegree, 'k-o', label='Hessenberg')
ax.loglog(degree, hstimeDegree[0]*(degree/degree[0])**2, 'k--', label='Quadratic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_hs_times_degree.pdf", format='pdf')

fig, ax = plt.subplots()
ax.loglog(degree, tdtimeDegree, 'k-o', label='Tridiagonal')
ax.loglog(degree, tdtimeDegree[0]*(degree/degree[0])**2, 'k--', label='Quadratic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_td_times_degree.pdf", format='pdf')

fig, ax = plt.subplots()
ax.loglog(size, getimeSize, 'k-o', label='General')
ax.loglog(size, getimeSize[0]*(size/size[0])**4, 'k--', label='Quartic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_ge_times_size.pdf", format='pdf')

fig, ax = plt.subplots()
ax.loglog(size, hstimeSize, 'k-o', label='Hessenberg')
ax.loglog(size, hstimeSize[0]*(size/size[0])**3, 'k--', label='Cubic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_hs_times_size.pdf", format='pdf')

fig, ax = plt.subplots()
ax.loglog(size, tdtimeSize, 'k-o', label='Tridiagonal')
ax.loglog(size, tdtimeSize[0]*(size/size[0])**2, 'k--', label='Quadratic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/cc_td_times_size.pdf", format='pdf')


