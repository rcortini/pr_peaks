#!/bin/bash

source ~/work/tools/my_env.sh

# returns the complete path of the BAM alignment file associated
# to the "sample_id" that user passes through the only argument
function chipseq_bam_location {
  sample_id=$1
  xavi_datadir='/mnt/xavi/data'
  d="$xavi_datadir/chipseq/samples/$sample_id/alignments"
  find $d -name "*.bam" | head -n 1
}

# check for proper invocation
if [ $# -lt 1 ]; then
  echo "remove_background.sh <input>" 1>&2
  exit 1
fi
input=$1

# build the file name for the input to use
if [ "$input" == "T0_roberto_input" ]; then
  input_bam=$(chipseq_bam_location $input)
else
  input_bam="/mnt/mbeato/rferrari/hg38/hg38_MMTV/comb_inputs/$input.bam"
fi

# check that the input file exists
if ! test -e $input_bam; then
  echo "Input file does not exist!" 1>&2
  exit 1
fi

# location of Zerone
zerone=$HOME/soft/zerone/zerone

# input file names
high_bam=$(chipseq_bam_location "gv_107_01_01_chipseq")
medium1_bam=$(chipseq_bam_location "gv_108_01_01_chipseq")
medium2_bam=$(chipseq_bam_location "gv_109_01_01_chipseq")
medium3_bam=$(chipseq_bam_location "gv_110_01_01_chipseq")
low_bam=$(chipseq_bam_location "gv_111_01_01_chipseq")

# output file names
data_dir=$HOME/work/CRG/projects/pr_peaks/data
output_dir=$data_dir/$input
mkdir -p $output_dir
high_out=$output_dir/high-zerone.out
medium1_out=$output_dir/medium1-zerone.out
medium2_out=$output_dir/medium2-zerone.out
medium3_out=$output_dir/medium3-zerone.out
low_out=$output_dir/low-zerone.out

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
