import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree1, dsTime, dsRad, pzTime, pzRad = [], [], [], [], []
degree2, amTime, amRad = [], [], []
with open('../results/outputSpoly.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree1.append(row[0].strip())
    dsTime.append(row[1].strip())
    dsRad.append(row[2].strip())
    pzTime.append(row[3].strip())
    pzRad.append(row[4].strip())

with open('../results/outputAMVW.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree2.append(row[0].strip())
    amTime.append(row[1].strip())
    amRad.append(row[2].strip())

fig, ax = plt.subplots()
ax.semilogy(degree1, dsTime, 'k-*', label='DSLMPEP Time')
ax.semilogy(degree1, pzTime, 'k-.', label='PZEROS Time')
ax.semilogy(degree2, amTime, 'k--', label='AMVW Time')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/spoly_times.pdf")

fig, ax = plt.subplots()
ax.semilogy(degree1, dsRad, 'k-*', label='DSLMPEP Radius')
ax.semilogy(degree1, pzRad, 'k-.', label='PZEROS Radius')
ax.semilogy(degree2, amRad, 'k--', label='AMVW Radius')
ax.set_ylabel('Forward Error')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/spoly_ferrs.pdf")
