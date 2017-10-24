import csv
import numpy as np
from matplotlib import pyplot as plt
with open('results/testConvergence.csv') as f:
  data =[]
  highlights= [['white', 'white', 'white', 'white'],\
              ['white', 'white', 'white', 'white'],\
              ['white', 'white', 'white', 'grey'],\
              ['white', 'white', 'white', 'grey'],\
              ['white', 'white', 'white', 'grey'],\
              ['white', 'white', 'grey', 'white'],\
              ['white', 'grey', 'grey', 'white'],\
              ['white', 'grey', 'grey', 'white'],\
              ['white', 'grey', 'white', 'white'],\
              ['white', 'white', 'white', 'white']]   
  reader = csv.reader(f)
  chosenTests = [3,9,10] #0 for row headers
  headers= next(reader)
  headers = [""]+[headers[i] for i in chosenTests]
  it = 1
  for row in reader:
    data.append(["Iter %s" % (it)]+[row[i].strip() for i in chosenTests])
    it = it+1

  nrows, ncols = len(data), len(headers)
  hcell, wcell = 0.175, 0.8 # tweak as per your requirements
  hpad, wpad = 0.5, 0.5    
  fig=plt.figure(figsize=(ncols*wcell+wpad, (nrows+1)*hcell+hpad))
  ax = fig.add_subplot(111)
  ax.axis('tight')
  ax.axis('off')
  the_table = ax.table(cellText=data,colLabels=headers,loc='center', cellColours = highlights)
  the_table.auto_set_font_size(False)
  the_table.set_fontsize(5.5)
  plt.savefig("diagrams/testConvergence.pdf", format='pdf')

