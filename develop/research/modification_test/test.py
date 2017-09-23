import csv
import numpy as np
from matplotlib import pyplot as plt

with open('results/results.csv') as f:
  reader = csv.reader(f)
  headers = next(reader)
  degree, dslm_time, dslm_berr, dslm1_time, dslm1_berr, dslm2_time, dslm2_berr, dsam_time, dsam_berr = [], [], [], [], [], [], [], [], []
  for row in reader:
    degree.append(float(row[0].strip()))
    dslm_time.append(float(row[1].strip()))
    dslm_berr.append(float(row[2].strip()))     
    dslm1_time.append(float(row[3].strip()))
    dslm1_berr.append(float(row[4].strip()))
    dslm2_time.append(float(row[5].strip()))
    dslm2_berr.append(float(row[6].strip()))
    dsam_time.append(float(row[7].strip()))
    dsam_berr.append(float(row[8].strip()))

  degree=np.array(degree)
  dslm_time=np.array(dslm_time)
  dslm_berr=np.array(dslm_berr)
  dslm1_time=np.array(dslm1_time)
  dslm1_berr=np.array(dslm1_berr)
  dslm2_time=np.array(dslm2_time)
  dslm2_berr=np.array(dslm2_berr)
  dsam_time=np.array(dsam_time)
  dsam_time=np.array(dsam_time)
  dsam_berr=np.array(dsam_berr)

  fig, axarr = plt.subplots(2, sharex=True)
  axarr[0].loglog(degree, dslm_time, 'k-o',  label='DSLM')
  axarr[0].loglog(degree, dslm1_time, 'k--', label='DSLM1')
  axarr[0].loglog(degree, dslm2_time, 'k--', label='DSLM2')
  axarr[0].loglog(degree, dsam_time, 'k-*', label='DSAM')
  axarr[0].set_ylabel('Seconds')
  axarr[0].set_xlabel('Degree')
  axarr[0].set_title('Times')
  legend = axarr[0].legend(loc=0, shadow=True)
  axarr[1].loglog(degree, dslm_berr, 'k-o', label='DSLM')
  axarr[1].loglog(degree, dslm1_berr, 'k--', label='DSLM1')
  axarr[1].loglog(degree, dslm2_berr, 'k--', label='DSLM2')
  axarr[1].loglog(degree, dsam_berr, 'k-*', label='DSAM')
  axarr[1].set_ylabel('Avg. of Max Backward Error')
  axarr[1].set_xlabel('Degree')
  axarr[1].set_title('Berr')
  plt.savefig("diagrams/test.pdf", format='pdf')

plt.show()

