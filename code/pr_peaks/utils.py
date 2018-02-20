from mybiotools import parse_simple_bed, chipseq_bam_location
import pysam

def load_hcp_peaks (peaks_id,xavi_datadir='/mnt/xavi') :
    datadir = '%s/projects/gvicent/2017-01-23_characterisation_prbs_r5020_titration/tables'%(xavi_datadir)
    datafile = '%s/genomic_coordinates_by_peak_population_%s.bed'%(datadir,peaks_id)
    return parse_simple_bed(datafile)

class Condition :
    def __init__(self,name,peak_code,concentration,sample_id) :
        self.name = name
        self.peak_code = peak_code
        self.concentration = concentration
        self.sample_id = sample_id
        # load the peaks
        self.peaks = load_hcp_peaks(self.peak_code)
        # init the BAM file
        self.bam_file = chipseq_bam_location(sample_id)
        # init the pysam parser
        self.bam = pysam.AlignmentFile(self.bam_file)
    def peak_counts(self,peak) :
        chromosome,start,end = peak
        chromosome = str(chromosome)
        # use the BigWig parser to get the stats of the peak
        return self.bam.count(chromosome,start,end)
    def __del__(self) :
        self.bam.close()

def average_peak_counts(peaks,condition) :
    npeaks = peaks.size
    pcounts = np.zeros(npeaks)
    for i,peak in enumerate(peaks) :
        pcounts[i] = condition.peak_counts(peak)
    return pcounts.mean()
