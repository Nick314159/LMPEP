import csv
import numpy as np
import sys
from matplotlib import pyplot as plt
methods = ["dslm", "dsam"]

for method in methods:
  with open('results/'+method+'_conv_test_results.txt') as f:
    reader = csv.reader(f)
    error= []
    for row in reader:
       error.append(float(row[0].strip()))   
    
    error=np.array(error)
    error = [elem for elem in error if elem > sys.float_info.epsilon]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(1, len(error)+1), error, 'k.', label='Error')
    plt.subplots_adjust(bottom=.1, left=.16)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Error')
    plt.title("Convergence Test "+method )
    legend = ax.legend(loc=0, shadow=True)
    plt.savefig("diagrams/"+method+"_conv_test_.pdf", format='pdf')

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



