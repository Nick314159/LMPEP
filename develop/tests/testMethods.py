from csv import reader
from matplotlib import pyplot as plt

with open('results/testMethods.csv') as f:
    table = reader(f)
    headers = next(table)
    degree, dslm_time, dslm_berr, dslm1_time, dslm1_berr, dsam_time, dsam_berr = [], [], [], [], [], [], []
    for row in table:
        degree.append(float(row[0].strip()))
        dslm_time.append(float(row[1].strip()))
        dslm_berr.append(float(row[2].strip()))     
        dslm1_time.append(float(row[3].strip()))
        dslm1_berr.append(float(row[4].strip()))
        dsam_time.append(float(row[5].strip()))
        dsam_berr.append(float(row[6].strip()))



    #fig, axarr = plt.subplots(2, sharex=True)
    fig, axarr = plt.subplots()
    axarr.loglog(degree, dslm_time, 'b-o',  label='Mod. Laguerre')
    axarr.loglog(degree, dslm1_time, 'k--', label='Laguerre ')
    axarr.loglog(degree, dsam_time, 'r-*', label='Aberth')
    axarr.set_ylabel('Elapsed Time (s)')
    axarr.set_xlabel('Degree')
    axarr.set_title('Time')
    legend = axarr.legend(loc=0, shadow=True)
    #axarr[1].loglog(degree, dslm_berr, 'b-o', label='Mod. Laguerre ')
    #axarr[1].loglog(degree, dslm1_berr, 'k--', label='Laguerre ')
    #axarr[1].loglog(degree, dsam_berr, 'r-*', label='Aberth')
    #axarr[1].set_ylabel('Max Error*')
    #axarr[1].set_xlabel('Degree')
    #axarr[1].set_title('Berr')
    #legend = axarr[1].legend(loc=0, shadow=True)
    #fig.subplots_adjust(hspace=.4)
    plt.savefig("diagrams/testMethods.pdf", format='pdf')

    plt.show()
