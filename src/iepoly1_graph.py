import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, size, nptimeDegree, nrtimeDegree, nptimeSize, nrtimeSize = [], [], [], [], [], []
#lberr, lferr, eaberr, eaferr =, [], [], [], []

with open('../results/outputIepoly1Degree.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(row[0].strip())
    nrtimeDegree.append(row[2].strip())
    nptimeDegree.append(row[3].strip())
with open('../results/outputIepoly1Size.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    size.append(row[1].strip())
    nrtimeSize.append(row[2].strip())
    nptimeSize.append(row[3].strip())
    
fig, ax = plt.subplots()
plt.title('Newton\'s Poylgon And Numerical Range, Degree Vs. Time Size 2')
ax.semilogy(degree, nptimeDegree, 'r-*', label='Newton\'s Polygon Time')
ax.semilogy(degree, nrtimeDegree, 'b-^', label='Numerical Range Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/iepoly_times_degree.pdf")

fig, ax = plt.subplots()
plt.title('Newton\'s Poylgon And Numerical Range, Size Vs. Time Degree 2')
ax.semilogy(degree, nptimeSize, 'r-*', label='Newton\'s Polygon Time')
ax.semilogy(degree, nrtimeSize, 'b-^', label='Numerical Range Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/iepoly_times_size.pdf")

