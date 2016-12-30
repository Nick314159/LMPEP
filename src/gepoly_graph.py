import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, size, ltimeDegree, eatimeDegree, ltimeSize, eatimeSize = [], [], [], [], [], []
#lberr, lferr, eaberr, eaferr =, [], [], [], []

with open('../results/outputGepolyDegree.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(row[0].strip())
    #size.append(row[1].strip())
    ltimeDegree.append(row[2].strip())
    #lberr.append(row[3].strip())
    #lferr.append(row[4].strip())
    eatimeDegree.append(row[5].strip())
    #eaberr.append(row[6].strip())
    #eaferr.append(row[7].strip())
with open('../results/outputGepolySize.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    #degree.append(row[0].strip())
    size.append(row[1].strip())
    ltimeSize.append(row[2].strip())
    #lberr.append(row[3].strip())
    #lferr.append(row[4].strip())
    eatimeSize.append(row[5].strip())
    #eaberr.append(row[6].strip())
    #eaferr.append(row[7].strip())

fig, ax = plt.subplots()
ax.semilogy(degree, ltimeDegree, 'r-*', label='Laguerre Time')
ax.semilogy(degree, eatimeDegree, 'b-^', label='Ehrlich-Aberth Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/gepoly_times_degree.pdf")

fig, ax = plt.subplots()
ax.semilogy(size, ltimeSize, 'r-*', label='Laguerre Time')
ax.semilogy(size, eatimeSize, 'b-^', label='Ehrlich-Aberth Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/gepoly_times_size.pdf")

