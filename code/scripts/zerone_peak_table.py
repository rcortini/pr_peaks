from mybiotools import log_message
import numpy as np
import zerone
import sys, os

def parse_simple_bed (fname) :
    simple_bed_dtype = np.dtype([
                                ('chr','S10'),
                                ('start',np.int64),
                                ('end',np.int64)
                               ])
    return np.genfromtxt(fname,dtype=simple_bed_dtype)

def load_hcp_peaks (peaks_id,xavi_datadir='/users/mbeato/public-docs') :
    datadir = '%s/projects/gvicent/2017-01-23_characterisation_prbs_r5020_titration/tables'%(xavi_datadir)
    datafile = '%s/genomic_coordinates_by_peak_population_%s.bed'%(datadir,peaks_id)
    return parse_simple_bed(datafile)

# check for proper invocation
if len(sys.argv) < 3 :
    print "zerone_peak_table <sample_id> <peaks_id>"
    sys.exit(1)

# get input from command line
sample_id = sys.argv[1]
peaks_id = sys.argv[2]

# prepare file names, parse zerone output file
log_message('zerone_peak_table','Parsing Zerone output file')
datadir = '%s/work/CRG/projects/pr_peaks/data/T0_roberto_input'%(os.getenv('HOME'))
zerone_fname = '%s/%s-zerone.out'%(datadir,sample_id)
experiment = zerone.ZeroneOutput(zerone_fname)
peaks = load_hcp_peaks(peaks_id)

# init table
table = np.zeros(peaks.size,dtype=np.int32)

# fill table
log_message('zerone_peak_table','Filling peak table')
for i,peak in enumerate(peaks) :
    peaks = experiment.find_peak(peak['chr'],peak['start'],peak['end'])
    enrichment = peaks['enrichment'].sum()/peaks.size
    table[i] = enrichment

# save output file
log_message('zerone_peak_table','Done. Saving output file')
out_fname = '%s/%s-%s.npy'%(datadir,sample_id,peaks_id)
np.save(out_fname,table)
