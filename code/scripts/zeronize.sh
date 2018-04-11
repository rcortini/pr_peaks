#!/bin/bash

source ~/work/tools/my_env.sh

# returns the complete path of the BAM alignment file associated
# to the "sample_id" that user passes through the only argument
function chipseq_bam_location {
  sample_id=$1
  if [[ `hostname` == *"ant-login"* ]]; then
    datadir='/users/mbeato/projects/data'
  else
    datadir='/mnt/mbeato/projects/data'
  fi
  d="$datadir/chipseq/samples/$sample_id/alignments"
  find $d -name "*.bam" | head -n 1
}

# check for proper invocation
if [ $# -lt 2 ]; then
  echo "zeronize.sh <input> <sample_id>" 1>&2
  exit 1
fi
input=$1
sample_id=$2

# build the file name for the input to use
if [ "$input" == "T0_roberto_input" ]; then
  input_bam=$(chipseq_bam_location $input)
else
  input_bam="/mnt/mbeato/rferrari/hg38/hg38_MMTV/comb_inputs/$input.bam"
fi

# check that the input file exists
if ! test -e $input_bam; then
  echo "Input file $input_bam does not exist!" 1>&2
  exit 1
fi

# location of Zerone
zerone=$HOME/soft/zerone/zerone

# input file names
bam=$(chipseq_bam_location $sample_id)

# output file names
data_dir=$HOME/work/CRG/projects/pr_peaks/data
output_dir=$data_dir/$input
if ! test -e $output_dir; then
  mkdir -p $output_dir
fi
out=$output_dir/$sample_id-zerone.out

# check that the output file exists: if so, exit without doing anything
if test -e $out; then
  echo "Output file $out already exists, aborting"
  exit 0
fi

# invoke Zerone
log_message "Zeroning $sample_id with $input"
$zerone -0 $input_bam --chip $bam > $out
log_message "Complete"
