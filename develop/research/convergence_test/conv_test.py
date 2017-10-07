import csv
import numpy as np
from matplotlib import pyplot as plt
methods = ["dslm", "dsam"]

for method in methods:
  with open('results/'+method+'_conv_test_results.txt') as f:
    reader = csv.reader(f)
    error= []
    for row in reader:
       error.append(float(row[0].strip()))   
    
    error=np.array(error)
    
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



