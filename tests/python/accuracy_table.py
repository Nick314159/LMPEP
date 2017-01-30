import csv
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pyplot
data = []
ferrdblm, ferrdlm, ferrsblm,ferrwslm= [],[],[],[]
ferrdbq, ferrdq, ferrsbq,ferrwsq= [],[],[],[]
with open('../results/outputAccuracy.csv') as f:
  reader = csv.reader(f)
  headers= next(reader)
  for row in reader:
    data.append([row[0].strip(), row[1].strip(), row[2].strip(), row[3].strip(), row[4].strip()])
    
fig=plt.figure()
ax = fig.add_subplot(111)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
colLabels=("Problem", "Max LM-BERR", "Avg LM-BERR", "Max QUAD-BERR", "Avg QUAD-BERR")
the_table = ax.table(cellText=data, colLabels=colLabels, loc='center')
plt.tight_layout()
plt.savefig("../results/accuracy_table.pdf", format='pdf')

with open('../results/outputAccuracyFerr.txt') as f:
  reader = csv.reader(f)
  headers= next(reader)
  i = 1
  for row in reader:
    if row==[]:
      i = i + 1
    else:
     if i == 1:
      ferrdblm.append(row[0].strip())
      ferrdbq.append(row[1].strip())
       
     if i == 2:
      ferrdlm.append(row[0].strip())
      ferrdq.append(row[1].strip())
       
     if i == 3:
      ferrsblm.append(row[0].strip())
      ferrsbq.append(row[1].strip())
       
     if i == 4:
      ferrwslm.append(row[0].strip())
      ferrwsq.append(row[1].strip())
   
   # dampled_beam, Dirac, speaker_box, and wiresaw2
fig, ax = plt.subplots()
ax.plot(ferrdblm, ferrdbq, 'ko', label='Dampled Beam')
ax.set_yscale('log')
ax.set_xscale('log')
lims = [
  np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
  np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_ylabel('QUAD')
ax.set_xlabel('LM')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/accuracy_ferr_db.pdf", format='pdf')

fig, ax = plt.subplots()
ax.plot(ferrdlm, ferrdq, 'ko', label='Dirac')
ax.set_yscale('log')
ax.set_xscale('log')
lims = [
  np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
  np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_ylabel('QUAD')
ax.set_xlabel('LM')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/accuracy_ferr_d.pdf", format='pdf')

fig, ax = plt.subplots()
ax.plot(ferrsblm, ferrsbq, 'ko', label='Speaker Box')
ax.set_yscale('log')
ax.set_xscale('log')
lims = [
  np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
  np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_ylabel('QUAD')
ax.set_xlabel('LM')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/accuracy_ferr_sb.pdf", format='pdf')

fig, ax = plt.subplots()
ax.plot(ferrwslm, ferrwsq, 'ko', label='Wiresaw 2')
ax.set_yscale('log')
ax.set_xscale('log')
lims = [
  np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
  np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_ylabel('QUAD')
ax.set_xlabel('LM')
legend = ax.legend(loc=0, shadow=True)
savefig("../results/accuracy_ferr_ws.pdf", format='pdf')

   

    

