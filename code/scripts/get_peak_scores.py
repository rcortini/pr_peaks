import numpy as np
import pr_peaks
import os
import mybiotools as mbt
from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.motifs.matrix import PositionWeightMatrix

# define a function that returns the maximum binding score detected on a list of sequences
def get_max_seq_score(genome,peaks,pssm) :
    max_scores = []
    for peak in peaks :
        seq = genome[peak['chr']].seq[peak['start']:peak['end']]
        try :
            scores_f = pssm.calculate(seq)
            scores_b = pssm.reverse_complement().calculate(seq)
            max_scores.append(max(scores_f.max(),scores_b.max()))
        except MemoryError :
            print seq, peak
    return np.array(max_scores)

# load the PR binding motif matrix so that BioPython understands it
M = np.genfromtxt(os.getenv('HOME')+'/work/data/motif231.motif',comments='>')
Mdict = {}
for i,letter in enumerate(['A','C','G','T']) :
    Mdict[letter] = M[:,i]
pwm = PositionWeightMatrix(IUPAC.unambiguous_dna,Mdict)
pssm = pwm.log_odds()
motif_length = len(pwm['A'])

# first I load the genome
hg19_genome_file = os.getenv('HOME') + '/work/data/GRCh37.fasta'
h19 = SeqIO.index (hg19_genome_file,'fasta',alphabet=IUPAC.unambiguous_dna)

# load the PR binding motif matrix so that BioPython understands it
M = np.genfromtxt(os.getenv('HOME')+'/work/data/motif231.motif',comments='>')
Mdict = {}
for i,letter in enumerate(['A','C','G','T']) :
    Mdict[letter] = M[:,i]
pwm = PositionWeightMatrix(IUPAC.unambiguous_dna,Mdict)
pssm = pwm.log_odds()
motif_length = len(pwm['A'])

# load the peak data
high       = pr_peaks.Condition('high'   ,'all_treated',0.05,'gv_107_01_01_chipseq')
medium1    = pr_peaks.Condition('medium1','4HCP'       ,0.10,'gv_108_01_01_chipseq')
medium2    = pr_peaks.Condition('medium2','3HCP'       ,0.50,'gv_109_01_01_chipseq')
medium3    = pr_peaks.Condition('medium3','3HCP'       ,1.00,'gv_110_01_01_chipseq')
low        = pr_peaks.Condition('low'    ,'1HCP'       ,10.0,'gv_111_01_01_chipseq')

# then calculate the peak scores
datadir = '../../data'
for condition in [high,medium2,low] :
    mbt.log_message('get_max_seq_score','peak_type = %s'%condition.peak_code)
    condition.max_peak_scores = get_max_seq_score(h19,condition.peaks,pssm)
    np.save('%s/%s-peak_scores.npy'%(datadir,condition.peak_code),condition.max_peak_scores)
