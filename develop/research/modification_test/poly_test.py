import csv
import numpy as np
from matplotlib import pyplot as plt

with open('results/poly_test_results.txt') as f:
  reader = csv.reader(f)
  temp = next(reader)
  while True: 
    test = temp[0].strip()
    dslm, dslm1, dslm2, dsam = [], [], [], []
    next(reader)
    for row in reader:
       if row[0].strip()=="DSLM1 Absolute Error:":
         break
       dslm.append(float(row[0].strip()))
    for row in reader:
       if row[0].strip()=="DSLM2 Absolute Error:":
         break
       dslm1.append(float(row[0].strip()))
    for row in reader:
       if row[0].strip()=="DSAM Absolute Error:":
         break
       dslm2.append(float(row[0].strip()))
    for row in reader:
       if row[0].strip().startswith("Deg"):
         temp = row
         break
       dsam.append(float(row[0].strip()))
    dslm=np.array(dslm)
    dslm1=np.array(dslm1)
    dslm2=np.array(dslm2)
    dsam=np.array(dsam)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(1, len(dslm1)+1), dslm1, 'w^', label='dslm1')
    ax.plot(range(1, len(dslm2)+1), dslm2, 'wo', label='dslm2')
    ax.plot(range(1, len(dslm)+1), dslm, 'k^', label='dslm')
    ax.plot(range(1, len(dsam)+1), dslm, 'ko', label='dsam')
    plt.subplots_adjust(bottom=.1, left=.16)
    ax.set_yscale('log')
    ax.set_xlabel('Eigenvalue Index')
    ax.set_ylabel('Relative Forward Error')
    plt.title(test)
    legend = ax.legend(loc=0, shadow=True)
    plt.savefig("diagrams/poly_test_"+test.replace(' ', '_')+".pdf", format='pdf')

plt.show()



