import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, size, ltimeDegree, eatimeDegree, ltimeSize, eatimeSize = [], [], [], [], [], []

with open('../results/outputGepolyDegree.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(row[0].strip())
    ltimeDegree.append(row[2].strip())
    eatimeDegree.append(row[5].strip())
with open('../results/outputGepolySize.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    size.append(row[1].strip())
    ltimeSize.append(row[2].strip())
    eatimeSize.append(row[5].strip())

fig, ax = plt.subplots()
ax.semilogy(degree, ltimeDegree, 'k-*', label='Laguerre Method')
ax.semilogy(degree, eatimeDegree, 'k--', label='Ehrlich-Aberth Method')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/gepoly_times_degree.pdf", format='pdf')

fig, ax = plt.subplots()
ax.semilogy(size, ltimeSize, 'k-*', label='Laguerre Method')
ax.semilogy(size, eatimeSize, 'k--', label='Ehrlich-Aberth Method')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/gepoly_times_size.pdf", format='pdf')

