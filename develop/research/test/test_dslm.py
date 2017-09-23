import csv
import numpy as np
from matplotlib import pyplot as plt

with open('results/results.csv') as f:
  reader = csv.reader(f)
  headers = next(reader)
  if headers[0].strip()=='DSLM':
    berr, degree, time = [], [], []
    for row in reader:
      degree.append(float(row[0].strip()))
      time.append(float(row[1].strip()))
      berr.append(float(row[2].strip()))

    degree=np.array(degree)
    time=np.array(time)
    berr=np.array(berr)
    fig, axarr = plt.subplots(2, sharex=True)
    axarr[0].loglog(degree, time, 'k-o',  label='Method Time')
    axarr[0].loglog(degree, time[0]*(degree/degree[0])**2, 'k--', label='O(Degree^2)')
    axarr[0].set_ylabel('Seconds')
    axarr[0].set_xlabel('Degree')
    axarr[0].set_title(headers[0])
    legend = axarr[0].legend(loc=0, shadow=True)
    axarr[1].loglog(degree, berr, 'k-o')
    axarr[1].set_ylabel('Avg. of Max Backward Error')
    axarr[1].set_xlabel('Degree')
    
      
  else:
    degree, time = [], []
    for row in reader:
      degree.append(float(row[0].strip()))
      time.append(float(row[1].strip()))

    degree=np.array(degree)
    time=np.array(time)
    fig, ax = plt.subplots()
    ax.loglog(degree, time, 'k-o',  label='Method Time')
    ax.loglog(degree, time[0]*(degree/degree[0]), 'k--', label='O(Degree)')
    ax.set_ylabel('Seconds')
    ax.set_xlabel('Degree')
    ax.set_title(headers[0])
    legend = ax.legend(loc=0, shadow=True)
    plt.savefig("diagrams/test_dslm.pdf", format='pdf')

plt.show()

