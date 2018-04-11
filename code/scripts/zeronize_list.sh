#!/bin/bash

root_dir=$(pwd)
pbs_in="zeronize.pbs.in"
while read line; do
  set -- $line
  input=$1
  sample_id=$2
  dir=cluster/$input/$sample_id
  if ! test -e $dir; then
    mkdir -p $dir
  fi
  pbs_out=$dir/zeronize.pbs
  cat $pbs_in |
    sed -e s,@INPUT@,$input,g |
    sed -e s,@SAMPLE_ID@,$sample_id,g |
  tee > $pbs_out
  cd $dir
    echo "qsub zeronize.pbs"
  cd $root_dir
done < zerone_list.txt
