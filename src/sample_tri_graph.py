import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
data = []
dgTime, qepTime = [], []
with open('../results/outputSampleTri.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    dgTime.append(row[1].strip())
    qepTime.append(row[2].strip())
    
fig, ax = plt.subplots()
ax.plot(dgTime, qepTime, 'ko', label='Times')
ax.set_xlabel('DGTLMPEP')
ax.set_ylabel('QEP3D')

lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
legend = ax.legend(loc=0, shadow=True)
plt.savefig("../results/sample_tri_graph.pdf", format='pdf')

