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
    degree1.append(float(row[0].strip()))
    dsTime.append(float(row[1].strip()))
    dsRad.append(float(row[2].strip()))
    pzTime.append(float(row[3].strip()))
    pzRad.append(float(row[4].strip()))

with open('../results/outputAMVW.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree2.append(float(row[0].strip()))
    amTime.append(float(row[1].strip()))
    amRad.append(float(row[2].strip()))

degree1 = np.array(degree1)
degree2 = np.array(degree2)
dsTime = np.array(dsTime)
pzTime = np.array(pzTime)
amTime = np.array(amTime)
dsRad = np.array(dsRad)
pzRad = np.array(pzRad)
amRad = np.array(amRad)

avg = (dsTime[0]+pzTime[0]+amTime[0])/3

fig, ax = plt.subplots()
ax.loglog(degree1, dsTime, 'k-*', label='DSLMPEP')
ax.loglog(degree1, pzTime, 'k-.', label='PZEROS')
ax.loglog(degree2, amTime, 'k--', label='AMVW')
ax.loglog(degree1, avg*(degree1/degree1[1])**2, ':', color='grey',
label='Quadratic Complexity')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/spoly_times.pdf", format='pdf')

fig, ax = plt.subplots()
ax.loglog(degree1, dsRad, 'k-*', label='DSLMPEP')
ax.loglog(degree1, pzRad, 'k-.', label='PZEROS')
ax.loglog(degree2, amRad, 'k--', label='AMVW')
ax.set_ylabel('Forward Error')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/spoly_ferrs.pdf", format='pdf')
