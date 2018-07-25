#!/bin/bash

source ~/work/tools/my_env.sh

#==================================================================================================
# Created on: 2017-01-23
# Usage: ./overlap_peaks_multiway_two_filters.sh
# Author: Javier Quilez (GitHub: jaquol)
# Modified by: Ruggero Cortini (GitHub: rcortini)
# Modification date: 2018-06-18
# Goal: calculate overlap between ChIP-seqs allowing q-value and enrichment filtering
#==================================================================================================

# gv_106_01_01_chipseq	T47D T0 PR
# gv_107_01_01_chipseq	T47D T30 PR 0.05nM
# gv_108_01_01_chipseq	T47D T30 PR 0.10nM
# gv_109_01_01_chipseq	T47D T30 PR 0.50nM
# gv_110_01_01_chipseq	T47D T30 PR 1nM
# gv_111_01_01_chipseq	T47D T30 PR 10nM



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

# variables
analysis=peak_analysis
samples="gv_106_01_01_chipseq gv_107_01_01_chipseq gv_108_01_01_chipseq gv_109_01_01_chipseq gv_110_01_01_chipseq gv_111_01_01_chipseq gv_066_01_01_chipseq rz_010_01_01_chipseq rz_011_01_01_chipseq rz_012_01_01_chipseq rz_013_01_01_chipseq"
min_qval=5 	# expressed as -log10(x)
min_enrichment=4
process=overlap_peaks
peak_caller=macs2
version=hg38_mmtv
peak_calling_mode=with_control
sequencing_type=single_end

# Paths
PROJECT=$HOME/work/CRG/projects/pr_peaks
DATA=$PROJECT/data/chipseq
ANALYSIS=$DATA/$analysis
mkdir -p $ANALYSIS


#==================================================================================================
# COMMANDS
#==================================================================================================

# filter peaks list by q-value and enrichment
for s in $samples; do
  ibed=$(find $DATA/samples/$s -name "*chipseq_peaks.narrowPeak")
  echo $ibed
  # ibed=$DATA/samples/$s/peaks/$peak_caller/$version/$peak_calling_mode/$sequencing_type/${s}_peaks.narrowPeak
  tbed=$ANALYSIS/$s
  log_message "Filtering $s"
  awk -v min_qval=$min_qval -v min_enrichment=$min_enrichment '($7 > min_enrichment) && ($9 > min_qval)' $ibed | grep -v "chrM\|chrY\|chrUn\|luciferase" > $tbed
done

# overlap with peaks with HOMER
otsv_venn=$ANALYSIS/overlap_peaks_minqval${min_qval}_minenrichment${min_enrichment}_venn.tsv
otsv_regions=$ANALYSIS/overlap_peaks_minqval${min_qval}_minenrichment${min_enrichment}_peaks_overlap.tsv
log_message "Merging peaks"
cd $ANALYSIS
mergePeaks -d given $samples -venn $otsv_venn > $otsv_regions
# rm -f $samples
log_message "Done"
cd
