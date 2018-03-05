import numpy as np
import matplotlib.pyplot as plt
import os

concentrations = np.array([0.05,0.10,0.50,1.00,10.0])

# load data
pr_peaks_root_dir = '%s/work/CRG/projects/pr_peaks'%(os.getenv('HOME'))
data_dir = '%s/data'%(pr_peaks_root_dir)
figures_dir = '%s/figures'%(pr_peaks_root_dir)
Hpeaks_count = np.load('%s/Hpeaks_count.npy'%(data_dir))
Mpeaks_count = np.load('%s/Mpeaks_count.npy'%(data_dir))
Lpeaks_average = np.load('%s/Lpeaks_average.npy'%(data_dir))

# plot now the number of reads as a function of the concentration for each peak class
fig = plt.figure(figsize=(6,4))
plt.semilogx(concentrations,Hpeaks_count.mean(axis=0)/Lpeaks_average,'o--',
            label='High')
plt.semilogx(concentrations,Mpeaks_count.mean(axis=0)/Lpeaks_average,'^--',
            label='Medium')
plt.xlabel('Concentration of hormone [nM]',fontsize=24)
plt.ylabel('H to L ratio',fontsize=24)
plt.legend(loc='upper right')
fig.tight_layout()
fig.savefig('%s/H_to_L.pdf'%(figures_dir))
plt.show()

fig = plt.figure(figsize=(6,4))
plt.semilogx(concentrations,Hpeaks_count.mean(axis=0),'o--',
            label='High')
plt.semilogx(concentrations,Mpeaks_count.mean(axis=0),'^--',
            label='Medium')
plt.xlabel('Concentration of hormone [nM]',fontsize=24)
plt.ylabel('Number of reads',fontsize=24)
plt.legend(loc='upper left')
fig.tight_layout()
fig.savefig('%s/peaks_average_count.pdf'%(figures_dir))
plt.show()
