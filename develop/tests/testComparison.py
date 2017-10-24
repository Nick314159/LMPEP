import csv
import numpy as np
from matplotlib import pyplot as plt

with open('results/testComparison1.csv') as f:
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
  axarr[0].loglog(degree, lm_time, 'b-o',  label='DSLM')
  axarr[0].loglog(degree, pz_time, 'k--', label='PZEROS')
  axarr[0].loglog(degree, a_time, 'r-*', label='AMVW')
  axarr[0].set_ylabel('Elapsed Time')
  axarr[0].set_xlabel('Degree')

  legend = axarr[0].legend(loc=0, shadow=True)
  axarr[1].loglog(degree, lm_berr, 'b-o', label='DSLM')
  axarr[1].loglog(degree, pz_berr, 'k--', label='PZEROS')
  axarr[1].loglog(degree, a_berr, 'r-*', label='AMVW')
  axarr[1].set_ylabel('Backward Error')
  axarr[1].set_xlabel('Degree')

  fig.subplots_adjust(hspace=.4)
  plt.savefig("diagrams/testComparison1.pdf", format='pdf')

with open('results/testComparison2.csv') as f:
  data =[]
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    data.append([row[0].strip(), row[1].strip(), row[2].strip(), row[3].strip()])

  nrows, ncols = len(data)+1, len(headers)
  hcell, wcell = 0.175, 0.8 # tweak as per your requirements
  hpad, wpad = 0.5, 0.5    
  fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
  ax = fig.add_subplot(111)
  ax.axis('tight')
  ax.axis('off')
  the_table = ax.table(cellText=data,colLabels=headers,loc='center')
  the_table.auto_set_font_size(False)
  the_table.set_fontsize(5.5)
  plt.savefig("diagrams/testComparison2.pdf", format='pdf')

