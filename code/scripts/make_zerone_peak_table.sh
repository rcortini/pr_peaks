#!/bin/bash

pbs_in="zerone_peak_table.pbs.in"
input="T0_roberto_input"
root_dir=$(pwd)
for sample_id in $(cat sample_ids.txt); do
  for peak_id in $(cat peak_ids.txt); do
    dir=cluster/$input/$sample_id/$peak_id
    if ! test -e $dir; then
      mkdir -p $dir
    fi
    pbs_out=$dir/zerone_peak_table.pbs
    cat $pbs_in |
      sed -e s,@SID@,$sample_id,g |
      sed -e s,@PID@,$peak_id,g |
    tee > $pbs_out
    cd $dir
      qsub zerone_peak_table.pbs
    cd $root_dir
  done
done
