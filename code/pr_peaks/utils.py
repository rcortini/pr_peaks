from mybiotools import parse_simple_bed

def load_hcp_peaks (peaks_id,xavi_datadir='/mnt/xavi') :
    datadir = '%s/projects/gvicent/2017-01-23_characterisation_prbs_r5020_titration/tables'%(xavi_datadir)
    datafile = '%s/genomic_coordinates_by_peak_population_%s.bed'%(datadir,peaks_id)
    return parse_simple_bed(datafile)
