import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, dsTime, dsRad, pzTime, pzRad = [], [], [], [], []
with open('results/output.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(row[0].strip())
    dsTime.append(row[1].strip())
    dsRad.append(row[2].strip())
    pzTime.append(row[3].strip())
    pzRad.append(row[4].strip())

fig, ax = plt.subplots()
#ax.title('Degree Vs. Time')
ax.semilogy(degree, dsTime, 'r-*', label='DSLMPEP Time')
ax.semilogy(degree, pzTime, 'b-^', label='PZEROS Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("times.pdf")

fig, ax = plt.subplots()
#ax.title('Radius Vs. Time')
ax.semilogy(degree, dsRad, 'r-*', label='DSLMPEP Radius')
ax.semilogy(degree, pzRad, 'b-^', label='PZEROS Radius')
ax.set_ylabel('Forward Error')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("forwardErrs.pdf")
