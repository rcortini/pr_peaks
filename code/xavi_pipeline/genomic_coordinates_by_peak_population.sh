#!/bin/bash

source ~/work/tools/my_env.sh

#==================================================================================================
# Created on: 2017-01-30
# Usage: ./genomic_coordinates_by_peak_population.sh
# Author: Javier Quilez (GitHub: jaquol)
# Modified by: Ruggero Cortini (GitHub: rcortini)
# Modification date: 2018-06-18
# Goal: generate lists of genomic coordinates for each of the specified population of peaks
#==================================================================================================

#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

# variables
analysis=peak_analysis
regions="1HCP 2HCP 3HCP 4HCP all_treated"
process=genomic_coordinates_by_peak_population
min_qval=10 	# expressed as -log10(x)
peak_caller=macs2
data_type=chipseq
version=hg38_mmtv
sequencing_type=single_end
peak_calling_mode=with_control
min_qval=5 	# expressed as -log10(x)
min_enrichment=4

# Paths
PROJECT=$HOME/work/CRG/projects/pr_peaks
DATA=$PROJECT/data
ANALYSIS=$DATA/$analysis
mkdir -p $ANALYSIS

peaks_overlap=$ANALYSIS/overlap_peaks_minqval${min_qval}_minenrichment${min_enrichment}_peaks_overlap.tsv

# hash array relating a certain tag to the peaks category
declare -A peaks_category_tag
peaks_category_tag["all_treated"]="gv_107_01_01_chipseq|gv_108_01_01_chipseq|gv_109_01_01_chipseq|gv_110_01_01_chipseq|gv_111_01_01_chipseq|gv_066_01_01_chipseq"
peaks_category_tag["4HCP"]="gv_108_01_01_chipseq|gv_109_01_01_chipseq|gv_110_01_01_chipseq|gv_111_01_01_chipseq|gv_066_01_01_chipseq"
peaks_category_tag["3HCP"]="gv_109_01_01_chipseq|gv_110_01_01_chipseq|gv_111_01_01_chipseq|gv_066_01_01_chipseq"
peaks_category_tag["2HCP"]="gv_110_01_01_chipseq|gv_111_01_01_chipseq|gv_066_01_01_chipseq"
peaks_category_tag["1HCP"]="gv_111_01_01_chipseq|gv_066_01_01_chipseq"

#==================================================================================================
# COMMANDS
#==================================================================================================

# prepare regions
for r in $regions; do

  log_message "Processing $r"
  obed=$ANALYSIS/genomic_coordinates_by_peak_population_$r.bed
  set_name=${peaks_category_tag[$r]}
  grep -v "#name" $peaks_overlap |awk -v var="$set_name" '$7 == var {OFS="\t"; print $2,$3-1,$4}' > $obed

done
