import csv
import numpy as np
from matplotlib import pyplot as plt

with open('results/testStart.csv') as f:
  reader = csv.reader(f)
  headers = next(reader)
  degree, cStartTime, bStartTime = [], [], []
  for row in reader:
    degree.append(float(row[0].strip()))
    cStartTime.append(float(row[1].strip()))
    bStartTime.append(float(row[2].strip()))     

  degree=np.array(degree)
  cStartTime=np.array(cStartTime)
  bStartTime=np.array(bStartTime)
  

  fig, axarr = plt.subplots()
  axarr.plot(degree, cStartTime, 'b-o',  label='Cameron Start Time')
  axarr.plot(degree, bStartTime, 'r-o', label='Bini Start Time')
  axarr.set_ylabel('Elapsed Time')
  axarr.set_xlabel('Degree')
  legend = axarr.legend(loc=0, shadow=True)

  fig.subplots_adjust(hspace=.4)
  plt.savefig("diagrams/testStart.pdf", format='pdf')

