import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
degree, dsTime, pzTime = [], [], []
with open('results/output.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    degree.append(row[0].strip())
    dsTime.append(row[1].strip())
    pzTime.append(row[2].strip())

fig, ax = plt.subplots()
ax.semilogy(degree, dsTime, 'r-')
ax.semilogy(degree, pzTime, 'b-')
#TODO legend
savefig("graph.pdf")

#image = pyplot.figure()
#pyplot.yscale('log')
#pyplot.ylabel('log(Time)')
#pyplot.xlabel('Degree')
#pyplot.plot(degree, dsTime, 'ro')
#pyplot.plot(degree, pzTime, 'bo')
#filename="graph.pdf"
#image.savefig(filename)
