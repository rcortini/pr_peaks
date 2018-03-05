import numpy as np
import mybiotools as mbt
import pr_peaks
import os

# load the peak data
high       = pr_peaks.Condition('high'   ,'all_treated',0.05,'gv_107_01_01_chipseq')
medium1    = pr_peaks.Condition('medium1','4HCP'       ,0.10,'gv_108_01_01_chipseq')
medium2    = pr_peaks.Condition('medium2','3HCP'       ,0.50,'gv_109_01_01_chipseq')
medium3    = pr_peaks.Condition('medium3','3HCP'       ,1.00,'gv_110_01_01_chipseq')
low        = pr_peaks.Condition('low'    ,'1HCP'       ,10.0,'gv_111_01_01_chipseq')
Hpeaks = high.peaks
Mpeaks = medium2.peaks
Lpeaks = low.peaks
conditions = [high,medium1,medium2,medium3,low]
nconditions = len(conditions)
nHpeaks = len(Hpeaks)
nMpeaks = len(Mpeaks)

# init the arrays
# Hpeaks_count = np.zeros((nHpeaks,nconditions))
Mpeaks_count = np.zeros((nMpeaks,nconditions))
averageL = np.zeros(nconditions)

# fill the arrays
for j,condition in enumerate(conditions) :
    mbt.log_message('peak_counts','condition = %s'%(condition.peak_code))
    # averageL[j] = pr_peaks.average_peak_counts(Lpeaks,condition)
    # for i,peak in enumerate(Hpeaks) :
        # Hpeaks_count[i,j] = condition.peak_counts(peak)
    for i,peak in enumerate(Mpeaks) :
        Mpeaks_count[i,j] = condition.peak_counts(peak)

# save data
datadir = '../../data'
# np.save('%s/Hpeaks_count.npy'%(datadir),Hpeaks_count)
np.save('%s/Mpeaks_count.npy'%(datadir),Mpeaks_count)
# np.save('%s/Lpeaks_average.npy'%(datadir),averageL)
