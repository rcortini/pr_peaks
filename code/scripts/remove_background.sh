#!/bin/bash

source ~/work/tools/my_env.sh

# returns the complete path of the BAM alignment file associated
# to the "sample_id" that user passes through the only argument
function chipseq_bam_location {
  sample_id=$1
  xavi_datadir='/mnt/xavi/data'
  d="$xavi_datadir/chipseq/samples/$sample_id/alignments"
  find $d -name "*.bam"
}

# location of Zerone
zerone=$HOME/soft/zerone/zerone

# input file names
input_bam=$(chipseq_bam_location "T0_roberto_input")
high_bam=$(chipseq_bam_location "gv_107_01_01_chipseq")
medium1_bam=$(chipseq_bam_location "gv_108_01_01_chipseq")
medium2_bam=$(chipseq_bam_location "gv_109_01_01_chipseq")
medium3_bam=$(chipseq_bam_location "gv_110_01_01_chipseq")
low_bam=$(chipseq_bam_location "gv_111_01_01_chipseq")

# output file names
data_dir=$HOME/work/CRG/projects/pr_peaks/data
high_out=$data_dir/high-zerone.out
medium1_out=$data_dir/medium1-zerone.out
medium2_out=$data_dir/medium2-zerone.out
medium3_out=$data_dir/medium3-zerone.out
low_out=$data_dir/low-zerone.out

# invoke Zerone
log_message "Zeroning HIGH"
$zerone -0 $input_bam --chip $high_bam > $high_out
log_message "Zeroning MEDIUM 1"
$zerone -0 $input_bam --chip $medium1_bam > $medium1_out
log_message "Zeroning MEDIUM 2"
$zerone -0 $input_bam --chip $medium2_bam > $medium2_out
log_message "Zeroning MEDIUM 3"
$zerone -0 $input_bam --chip $medium3_bam > $medium3_out
log_message "Zeroning LOW"
$zerone -0 $input_bam --chip $low_bam > $low_out
