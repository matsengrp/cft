#!/usr/bin/env bash

from="$PARTIS/test/reference-results"
to="sampledata/QW333.043-Vk"

#cp $from/partition-new-simu-cluster-annotations.csv $to/Hs-LN2-5RACE-IgG-new-cluster-annotations.csv
#cp $from/partition-new-simu.csv $to/Hs-LN2-5RACE-IgG-new.csv


#bin/process_partis.py \
   #--annotations $from/partition-ref-simu-cluster-annotations.csv \
   #--partition $from/partition-ref-simu.csv \
   #--partis_log $from/test.log \
   #--cluster_base cluster \
   #--output_dir output_post_partis


# This works running manually
#p=/fh/fast/matsen_e/processed-data/partis/kate-qrs-2016-09-09/seeds/QA255.006-Vh
#pref=Hs-LN2-5RACE-IgG

#echo
  #"bin/process_partis.py \
     #--annotations $p/$pref-cluster-annotations.csv \
     #--partition $p/$pref.csv \
     #--partis_log $p/$pref.log \
     #--cluster_base cluster \
     #--output_dir output_post_partis"

