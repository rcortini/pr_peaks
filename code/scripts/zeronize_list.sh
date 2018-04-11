#!/bin/bash

root_dir=$(pwd)
while read line; do
  set -- $line
  input=$1
  sample_id=$2
  dir=cluster/$input/$sample_id
  pbs_out=$dir/zeronize.pbs
  cat $pbs_out |\
    sed -e s,@INPUT@,$input,g |\
    sed -e s,@SAMPLE_ID@,$sample_id,g |\
    cd $dir
      echo "qsub $pbs_out"
    cd $root_dir
done < zerone_list.txt
