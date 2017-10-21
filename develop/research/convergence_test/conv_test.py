import csv
import numpy as np
from matplotlib import pyplot as plt
with open('results/dslm_conv_test_results.txt') as f:
  reader = csv.reader(f)
  data, error, cols= [], [], []
  rows = [i for i in range(61)]
  for row in reader:
     if row[0].strip().startswith("Deg"):
       cols.append(str(row[0].strip()))
       if(len(error)>0):
         data.append(error)
       error = []  
     else:
       error.append(float(row[1].strip()))
 
  data.append(error)
  data = zip(*data[::-1])
  nrows, ncols = len(data)+1, len(cols)
  hcell, wcell = 0.175, 0.75 # tweak as per your requirements
  hpad, wpad = 0.5, 0.5    
  fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
  ax = fig.add_subplot(111)
  ax.axis('tight')
  ax.axis('off')
  the_table = ax.table(cellText=data,colLabels=cols,loc='center', rowLabels=rows)
  the_table.auto_set_font_size(False)
  the_table.set_fontsize(5.5)
  plt.savefig("diagrams/dslm_conv_test_.pdf", format='pdf')

