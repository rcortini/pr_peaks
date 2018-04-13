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

def zerone_outfile (sample_id) :
    datadir = '/users/project/prj005866_no_backup/pr_peaks/data/T0_roberto_input'
    return '%s/%s-zerone.out'%(datadir,sample_id)

# check for proper invocation
if len(sys.argv) < 3 :
    print "zerone_peak_table <sample_id> <peaks_id>"
    sys.exit(1)

# get input from command line
sample_id = sys.argv[1]
peaks_id = sys.argv[2]

zerone_fname = zerone_outfile(sample_id)
print os.path.exists(zerone_fname)
peaks = load_hcp_peaks(peaks_id)
