import csv
import numpy as np
from matplotlib import pyplot as plt
with open('results/testConvergence.csv') as f:
  data =[]
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    data.append([row[0].strip(),row[1].strip(),row[2].strip()])

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
  plt.savefig("diagrams/testConvergence.pdf", format='pdf')

