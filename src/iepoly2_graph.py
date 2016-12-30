import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
test, npRealEstimate, npImaginaryEstimate, nrRealEstimate, nrImaginaryEstimate, npReal, npImaginary, nrReal, nrImaginary = [], [], [], [], [], [], [], [], []


with open('../results/outputIepoly2.txt') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(row[0].strip())
    nrtimeDegree.append(row[2].strip())
    nptimeDegree.append(row[3].strip())

fig, ax = plt.subplots()
ax.semilogy(degree, nptimeDegree, 'r-*', label='Newtons Polygon Time')
ax.semilogy(degree, nrtimeDegree, 'b-^', label='Numerical Range Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/gepoly_times_degree.pdf")

fig, ax = plt.subplots()
ax.semilogy(degree, nptimeDegree, 'r-*', label='Newtons Polygon Time')
ax.semilogy(degree, nrtimeDegree, 'b-^', label='Numerical Range Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/gepoly_times_size.pdf")

