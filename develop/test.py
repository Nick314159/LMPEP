import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, time = [], []

with open('results.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(float(row[0].strip()))
    time.append(float(row[1].strip()))

degree=np.array(degree)
time=np.array(time)
fig, ax = plt.subplots()
ax.loglog(degree, time, 'k-o',  label='Method Time')
ax.loglog(degree, time[0]*(degree/degree[0]), 'k--', label='O(n)')
ax.set_ylabel('Seconds')
ax.set_xlabel('Degree')
legend = ax.legend(loc=0, shadow=True)
savefig("results.pdf", format='pdf')

