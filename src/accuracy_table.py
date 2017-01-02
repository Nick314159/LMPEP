import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
problem, berr, ferr, data = [], [], [], []

with open('../results/outputAccuracy.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    problem.append(row[0].strip())
    berr.append(row[1].strip())
    ferr.append(row[2].strip())
    
fig=plt.figure()
ax = fig.add_subplot(111)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
colLabels=("Problem", "BERR", "FERR")
data = [problem, berr, ferr]
the_table = ax.table(cellText=data,
          colLabels=colLabels,
          loc='center')
plt.savefig("accuracy_table.png")

