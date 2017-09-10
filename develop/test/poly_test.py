import csv
import numpy as np
from matplotlib import pyplot as plt

with open('poly_test_results.txt') as f:
  reader = csv.reader(f)
  temp = next(reader)
  while True: 
    test = temp[0].strip()
    lmpep, pzeros, amvw = [], [], []
    next(reader)
    for row in reader:
       if row[0].strip()=="POLZEROS Absolute Error:":
         break
       lmpep.append(float(row[0].strip()))
    for row in reader:
       if row[0].strip()=="AMVW Absolute Error:":
         break
       pzeros.append(float(row[0].strip()))
    for row in reader:
       if row[0].strip().startswith("Deg"):
         temp = row
         break
       amvw.append(float(row[0].strip()))
    lmpep=np.array(lmpep)
    pzeros=np.array(pzeros)
    amvw=np.array(amvw)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(1, len(lmpep)+1), lmpep, 'k*', label='LMPEP')
    ax.plot(range(1, len(pzeros)+1), pzeros, 'w^', label='PZEROS')
    ax.plot(range(1, len(amvw)+1), amvw, 'wo', label='AMVW')
    plt.subplots_adjust(bottom=.1, left=.16)
    ax.set_yscale('log')
    ax.set_xlabel('Eigenvalue Index')
    ax.set_ylabel('Relative Forward Error')
    plt.title(test)
    legend = ax.legend(loc=0, shadow=True)
    plt.savefig("poly_test_"+test.replace(' ', '_')+".pdf", format='pdf')

plt.show()



