import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
size, triTime, qepTime = [], [], []
with open('../results/outputSpoly.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    size.append(row[0].strip())
    triTime.append(row[1].strip())
    qepTime.append(row[2].strip())

fig, ax = plt.subplots()
ax.loglog(size, triTime, 'k-*', label='DGTLMPEP')
ax.loglog(size, qepTime, 'k-.', label='QEP3D')
ax.set_ylabel('Seconds')
ax.set_xlabel('Size')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/tri_times.pdf", format='pdf')
