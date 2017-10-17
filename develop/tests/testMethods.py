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



    fig, axarr = plt.subplots(2, sharex=True)
    axarr[0].loglog(degree, dslm_time, 'b-o',  label='Mod. Laguerre')
    axarr[0].loglog(degree, dslm1_time, 'k--', label='Laguerre ')
    axarr[0].loglog(degree, dsam_time, 'r-*', label='Aberth')
    axarr[0].set_ylabel('Seconds')
    axarr[0].set_xlabel('Degree')
    axarr[0].set_title('Elapsed Time')
    legend = axarr[0].legend(loc=0, shadow=True)
    axarr[1].loglog(degree, dslm_berr, 'b-o', label='Mod. Laguerre ')
    axarr[1].loglog(degree, dslm1_berr, 'k--', label='Laguerre ')
    axarr[1].loglog(degree, dsam_berr, 'r-*', label='Aberth')
    axarr[1].set_ylabel('Avg. of Max Backward Error')
    axarr[1].set_xlabel('Degree')
    axarr[1].set_title('Berr')
    legend = axarr[1].legend(loc=0, shadow=True)
    fig.subplots_adjust(hspace=.5)
    plt.savefig("diagrams/testMethods.pdf", format='pdf')

    plt.show()
