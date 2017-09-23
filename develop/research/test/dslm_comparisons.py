import csv
import numpy as np
from matplotlib import pyplot as plt

with open('results/results.csv') as f:
  reader = csv.reader(f)
  headers = next(reader)
  degree, lm_time, lm_berr, pz_time, pz_berr, a_time, a_berr = [], [], [], [], [], [], []
  for row in reader:
    degree.append(float(row[0].strip()))
    lm_time.append(float(row[1].strip()))
    lm_berr.append(float(row[2].strip()))     
    pz_time.append(float(row[3].strip()))
    pz_berr.append(float(row[4].strip()))
    a_time.append(float(row[5].strip()))
    a_berr.append(float(row[6].strip()))

  degree=np.array(degree)
  lm_time=np.array(lm_time)
  lm_berr=np.array(lm_berr)
  pz_time=np.array(pz_time)
  pz_berr=np.array(pz_berr)
  a_time=np.array(a_time)
  a_berr=np.array(a_berr)

  fig, axarr = plt.subplots(2, sharex=True)
  axarr[0].loglog(degree, lm_time, 'k-o',  label='LMPEP')
  axarr[0].loglog(degree, pz_time, 'k--', label='Pzeroes')
  axarr[0].loglog(degree, a_time, 'k-*', label='AMVW')
  axarr[0].set_ylabel('Seconds')
  axarr[0].set_xlabel('Degree')
  axarr[0].set_title('Times')
  legend = axarr[0].legend(loc=0, shadow=True)
  axarr[1].loglog(degree, lm_berr, 'k-o', label='LMPEP')
  axarr[1].loglog(degree, pz_berr, 'k--', label='Pzeroes')
  axarr[1].loglog(degree, a_berr, 'k-*', label='AMVW')
  axarr[1].set_ylabel('Avg. of Max Backward Error')
  axarr[1].set_xlabel('Degree')
  axarr[1].set_title('Berr')
  plt.savefig("diagrams/dslm_comparisions".pdf", format='pdf')

plt.show()

