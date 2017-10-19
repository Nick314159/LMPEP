import csv
import numpy as np
import sys
import random
from matplotlib import pyplot as plt
from matplotlib import colors 
methods = ["dslm"]

with open('results/dslm_conv_test_results.txt') as f:
  reader = csv.reader(f)
  error= []
  i = 0;
  temp = []
  for row in reader:
     if len(row)>0:
       temp.append(float(row[0].strip()))   
     else:
       i+=1
       error.append(temp)
       temp=[]
    
  fig = plt.figure()
  ax = fig.add_subplot(111)
  for j in range(i):
    temp=np.array(error[j])
    temp = [elem for elem in temp if elem > sys.float_info.epsilon]
    ax.scatter(range(1, len(temp)+1), temp, marker=".", color="#%06x" % random.randint(0, 0xFFFFFF))
  plt.subplots_adjust(bottom=.1, left=.16)
  ax.set_yscale('log')
  #ax.set_xscale('log')
  ax.set_xlabel('Iteration')
  ax.set_ylabel('Error')
  plt.title("Convergence Test DSLM" )
  legend = ax.legend(loc=0, shadow=True)
  plt.savefig("diagrams/dslm_conv_test_.pdf", format='pdf')

    #x = error[1:len(error)-1]
    #y = error[2:len(error)]
    #m,b = np.polyfit(x, y, 1)
    #l = [m * i + b for i in x]
  
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(x, y, 'k.', label='Error')
      
    #x_new = np.linspace(x[0], x[-1], num=len(x)*10)
    ##coeff = np.polyfit(x, y, 1)
    ##fit = np.poly1d(coeff)
    ##plt.plot(x_new, fit(x_new))
    #ax.plot(x, l, 'k', label="y={}x+{}".format(m,b))  

    #plt.subplots_adjust(bottom=.1, left=.16)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #ax.set_xlabel('Error at Nth')
    #ax.set_ylabel('Error at N+1st')
    #plt.title("Convergence Test "+method )
    #legend = ax.legend(loc=0, shadow=True)
    #plt.savefig("diagrams/"+method+"_filtered_conv_test_.pdf", format='pdf')



