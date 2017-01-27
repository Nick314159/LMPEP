import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
data = []
with open('../results/outputTri.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    data.append([row[0].strip(),row[1].strip(),row[2].strip(),row[3].strip(),row[4].strip()])
   
fig=plt.figure()
ax = fig.add_subplot(111)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
colLabels=("Problem", "DGTLMPEP TIME", "DGTLMPEP AVG. FERR", "EIGEN TIME", "EIGEN AVG. FERR")
the_table = ax.table(cellText=data, colLabels=colLabels, loc='center')
plt.tight_layout()
savefig("../results/tri_times_table.pdf", format='pdf')
