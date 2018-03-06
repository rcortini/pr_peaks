import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# directory names
pr_peaks_root_dir = '%s/work/CRG/projects/pr_peaks'%(os.getenv('HOME'))
data_dir = '%s/data'%(pr_peaks_root_dir)
figures_dir = '%s/figures'%(pr_peaks_root_dir)

# prepare labels
peaks = {'all_treated':'High',
         '3HCP':'Medium',
         '1HCP':'Low'}

# prepare figure
fig = plt.figure(figsize=(6,4))
x = np.arange(1,16,0.1)

# load data and calculate Gaussian kernels
for peak_name,label in peaks.iteritems() :
    data = np.load('%s/%s-peak_scores.npy'%(data_dir,peak_name))
    mask = ~np.isnan(data)
    kernel = gaussian_kde(data[mask])
    plt.plot(x,kernel(x),linewidth=3,label=label)

# finalize
xticks = range(2,18,2)
plt.xticks(xticks)
plt.xlabel('Site score',fontsize=24)
plt.ylabel('Distribution',fontsize=24)
plt.legend(loc='upper right')
fig.tight_layout()
fig.savefig('%s/peak_scores.pdf'%(figures_dir))
plt.show()
